import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import io
import re

# ==========================================
# 0. 网页配置与参数设置
# ==========================================
st.set_page_config(page_title="FeliEase 药物筛选分析平台 (激活剂模式)", layout="wide")

st.title("📈 FeliEase 高通量药物筛选分析平台 - 激活剂/升高信号筛选")
st.markdown("""
此工具用于自动化分析药物筛选数据，**寻找引起信号升高的化合物**。
- **自动定位板位**：彻底无视空行，自动巡航识别每板位置。
- **智能清洗**：自动剔除内参异常值和毒性假阳性。
- **化合物映射**：请上传 MCE 的 **Detailed Information 列表**，自动为您匹配药物真实姓名！
""")

# 侧边栏：参数配置
with st.sidebar:
    st.header("⚙️ 参数配置")
    st.success("🤖 自动找板引擎已激活\n无需手动设置空行间距！")
    PLATE_HEIGHT = st.number_input("每板行数 (Data Rows)", value=8, help="例如 B1-M8 是8行")
    
    st.divider()
    THRESHOLD_HIT = st.slider("Hit 判定阈值 (Ratio ≥ 该值)", 1.0, 10.0, 2.0, 0.5)
    LIMIT_HIGH_SIGNAL = st.number_input("高信号异常剔除倍数 (防沉淀/杂质)", value=15.0)
    
    st.divider()
    st.info("默认配置：\n内参: 第2列 (前5个)\n毒性: 第13列 (前5个)\n药物: 第3-12列")

COL_IDX_PLATE_ID = 0
COL_IDX_CONTROL = 1
COL_IDX_TOX = 12
COLS_IDX_DRUG = range(2, 12)
ROWS_CONTROL_RELATIVE = range(0, 5)
ROWS_TOX_RELATIVE = range(0, 5)

# ==========================================
# 1. 核心计算与辅助函数
# ==========================================
def calculate_clean_mean(values):
    data = pd.to_numeric(values, errors='coerce')
    data = data[~np.isnan(data)]
    if len(data) < 3: return np.mean(data) if len(data)>0 else 0
    mean_val = np.mean(data)
    std_val = np.std(data, ddof=1)
    clean = data[(data >= mean_val - 2*std_val) & (data <= mean_val + 2*std_val)]
    return np.mean(clean)

def get_col_letter(col_idx):
    if 0 <= col_idx <= 25: return chr(65 + col_idx)
    return f"C{col_idx}"

def format_well(w):
    """ 暴力清洗孔位：切除所有空格、Tab、换行符及非字母数字字符 """
    w = str(w)
    w = re.sub(r'[^a-zA-Z0-9]', '', w)
    if len(w) >= 2 and w[0].isalpha() and w[1:].isdigit():
        return f"{w[0].upper()}{int(w[1:]):02d}"
    return w

def load_mce_library(file_obj, filename):
    """ 超级健壮的化合物库读取器，完美处理回车、空行和乱码 """
    try:
        # 1. 安全读取（防回车截断）
        if filename.lower().endswith('.csv'):
            try:
                df = pd.read_csv(file_obj, header=None, dtype=str, encoding='utf-8-sig', on_bad_lines='skip')
            except:
                file_obj.seek(0)
                df = pd.read_csv(file_obj, header=None, dtype=str, encoding='gbk', on_bad_lines='skip')
        else:
            df = pd.read_excel(file_obj, header=None, dtype=str)
            
        # 2. 智能搜索表头
        header_idx = -1
        for i in range(min(50, len(df))):
            row_str = " ".join(df.iloc[i].fillna("").astype(str)).lower()
            if "plate" in row_str and "well" in row_str:
                header_idx = i
                break
                
        if header_idx != -1:
            raw_headers = df.iloc[header_idx].fillna("").astype(str)
            cleaned_headers = [re.sub(r'[\x00-\x1F\x7F\ufeff]', '', h).strip() for h in raw_headers]
            df.columns = cleaned_headers
            df = df.iloc[header_idx+1:].reset_index(drop=True)
        else:
            st.warning("⚠️ 未找到 Plate 和 Well 列！您可能上传了矩阵状的 'Layout' 表，请上传长列表状的 'Detailed Information' 表！")
            return None
            
        # 3. 规范化列名
        col_map = {}
        for c in df.columns:
            cl = c.lower().strip()
            if cl in ['plate', 'rack', 'platenum']: col_map[c] = 'Plate'
            elif cl in ['well', 'position']: col_map[c] = 'Well'
            elif 'product name' in cl or 'compound name' in cl or cl == 'name' or 'productname' in cl: col_map[c] = 'Product Name'
            elif 'catalog' in cl or 'item' in cl: col_map[c] = 'Catalog Number'
            elif 'target' in cl: col_map[c] = 'Target'
        df.rename(columns=col_map, inplace=True)
        
        # 4. 彻底清洗垃圾数据
        if 'Well' in df.columns:
            df['Well'] = df['Well'].apply(lambda x: re.sub(r'[^a-zA-Z0-9]', '', str(x)))
            valid_well_mask = df['Well'].str.match(r'^[A-Za-z]\d{1,2}$')
            df = df[valid_well_mask].copy()
            
        if 'Plate' in df.columns:
            unique_plates = df['Plate'].dropna().unique()
            unique_plates = [p for p in unique_plates if str(p).strip() != '']
            plate_mapping = {p: str(i+1) for i, p in enumerate(unique_plates)}
            df['Physical_Plate'] = df['Plate'].map(plate_mapping)
            
        if 'Well' in df.columns:
            df['Physical_Well'] = df['Well'].apply(format_well)
            
        return df
    except Exception as e:
        st.error(f"解析化合物库失败: {e}")
        return None

# ==========================================
# 2. 主逻辑
# ==========================================
col_upload1, col_upload2 = st.columns(2)
with col_upload1:
    uploaded_file = st.file_uploader("1️⃣ 必填：请上传 CSV 筛药结果", type=["csv"])
with col_upload2:
    uploaded_lib = st.file_uploader("2️⃣ 可选：上传 MCE Detailed Information 表", type=["csv", "xlsx"])

if uploaded_file is not None:
    try:
        df = pd.read_csv(uploaded_file, header=None)
        
        all_results = []
        total_tox_drop = 0
        total_high_drop = 0
        global_drug_count = 0
        plate_count = 1  
        
        progress_bar = st.progress(0)
        
        # ---------------------------------------------
        # 【超级找板引擎】自动无视所有不规则空行！
        # ---------------------------------------------
        i = 0
        while i < len(df):
            val = str(df.iloc[i, COL_IDX_PLATE_ID]).strip()
            
            if val and val.lower() != 'nan':
                if i + int(PLATE_HEIGHT) <= len(df):
                    plate_df = df.iloc[i : i + int(PLATE_HEIGHT), :].reset_index(drop=True)
                    plate_id = val
                    physical_plate = str(plate_count) 
                    
                    try:
                        curr_ctrls = plate_df.iloc[ROWS_CONTROL_RELATIVE, COL_IDX_CONTROL].values
                        plate_ctrl_mean = calculate_clean_mean(curr_ctrls)
                        
                        curr_toxs = plate_df.iloc[ROWS_TOX_RELATIVE, COL_IDX_TOX].values
                        plate_tox_threshold = calculate_clean_mean(curr_toxs)
                    except:
                        plate_ctrl_mean = 0
                        
                    if plate_ctrl_mean > 0:
                        drug_block = plate_df.iloc[:, COLS_IDX_DRUG]
                        for (rel_r, abs_c), cell_val in drug_block.stack().items():
                            global_drug_count += 1
                            try: 
                                cell_val = float(cell_val)
                            except: 
                                continue
                                
                            if np.isnan(cell_val): continue

                            if cell_val < plate_tox_threshold:
                                total_tox_drop += 1
                                continue
                            if cell_val > (plate_ctrl_mean * LIMIT_HIGH_SIGNAL):
                                total_high_drop += 1
                                continue
                            
                            ratio = cell_val / plate_ctrl_mean
                            is_hit = "Yes" if ratio >= THRESHOLD_HIT else "No"
                            
                            excel_row_num = i + rel_r + 1
                            excel_col_char = get_col_letter(abs_c)
                            
                            plate_row_char = chr(65 + rel_r)
                            physical_well = f"{plate_row_char}{abs_c:02d}"
                            
                            all_results.append({
                                "Global_ID": global_drug_count,
                                "Plate_ID": plate_id,            
                                "Physical_Plate": physical_plate, 
                                "Physical_Well": physical_well,  
                                "Coordinate": f"{excel_col_char}{excel_row_num}",
                                "Ratio": ratio,
                                "Raw_RFU": cell_val,
                                "Is_Hit": is_hit
                            })
                    
                    plate_count += 1
                    i += int(PLATE_HEIGHT) 
                    progress_bar.progress(min(i / len(df), 1.0))
                else:
                    break 
            else:
                i += 1 

        progress_bar.progress(1.0)

        # ==========================
        # 合并化合物信息库
        # ==========================
        if all_results:
            df_res = pd.DataFrame(all_results)
            
            if uploaded_lib is not None:
                st.info("🔄 正在智能解析化合物库...")
                df_lib = load_mce_library(uploaded_lib, uploaded_lib.name)
                
                if df_lib is not None:
                    if "Physical_Plate" not in df_lib.columns or "Physical_Well" not in df_lib.columns:
                        st.warning("⚠️ 文件中未能成功生成 Plate 和 Well 对照列。")
                    else:
                        df_res['Physical_Plate'] = df_res['Physical_Plate'].astype(str)
                        df_res['Physical_Well'] = df_res['Physical_Well'].astype(str)
                        df_lib['Physical_Plate'] = df_lib['Physical_Plate'].astype(str)
                        df_lib['Physical_Well'] = df_lib['Physical_Well'].astype(str)
                        
                        df_res = pd.merge(df_res, df_lib, on=["Physical_Plate", "Physical_Well"], how="left")
                        
                        # 检查是否真有匹配上的
                        if 'Product Name' in df_res.columns and df_res['Product Name'].notna().sum() > 0:
                            st.success("✅ 化合物名称完美匹配！")
                        else:
                            st.warning("⚠️ 化合物库已读取，但没有任何孔位匹配成功。请检查筛药结果和 MCE 库是不是同一批次！")

            # --------------------------
            # 展示与绘图
            # --------------------------
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("总扫描孔数", global_drug_count)
            col2.metric("毒性/假阳性剔除", total_tox_drop)
            col3.metric("异常极高信号剔除", total_high_drop)
            col4.metric("有效数据量", len(df_res))
            
            st.divider()
            
            st.subheader("📊 筛选结果可视化")
            fig, ax = plt.subplots(figsize=(12, 6))
            ax.scatter(df_res["Global_ID"], df_res["Ratio"], c='#5cb85c', s=15, alpha=0.6, label='Candidates')
            
            hits_df = df_res[df_res["Is_Hit"] == "Yes"]
            if not hits_df.empty:
                ax.scatter(hits_df["Global_ID"], hits_df["Ratio"], c='red', s=40, label=f'Hits (≥ {THRESHOLD_HIT})')
                
            ax.axhline(1.0, color='black', linestyle='--', label='Plate Control (1.0)')
            ax.axhline(THRESHOLD_HIT, color='blue', linestyle=':', label=f'Hit Threshold ({THRESHOLD_HIT})')
            ax.set_xlabel('Global Drug Sequence')
            ax.set_ylabel('Normalized Ratio')
            
            y_max = df_res["Ratio"].max()
            ax.set_ylim(bottom=-0.1, top=max(y_max * 1.1, THRESHOLD_HIT * 1.5))
            ax.legend()
            ax.grid(True, linestyle=':', alpha=0.4)
            st.pyplot(fig)
            
            st.divider()
            st.subheader("🎉 强效升高信号药 (Hits) 列表")
            
            if not hits_df.empty:
                base_cols = ["Plate_ID", "Physical_Well", "Ratio", "Raw_RFU"] 
                lib_cols = [c for c in ["Product Name", "Catalog Number", "Target", "Pathway"] if c in df_res.columns]
                display_cols = lib_cols + base_cols
                
                hits_display = hits_df.sort_values(by="Ratio", ascending=False)
                final_cols = [c for c in display_cols if c in hits_display.columns]
                
                st.dataframe(hits_display[final_cols].head(15))
                
                buffer = io.BytesIO()
                with pd.ExcelWriter(buffer, engine='xlsxwriter') as writer:
                    export_cols = final_cols + [c for c in hits_display.columns if c not in final_cols and c not in ["Physical_Plate", "Is_Hit", "Coordinate"]]
                    hits_display[export_cols].to_excel(writer, index=False, sheet_name='Hits')
                
                st.download_button(
                    label="📥 下载带有化合物详细信息的 Hits (Excel)",
                    data=buffer.getvalue(),
                    file_name="Enhancer_Hits_with_Info.xlsx",
                    mime="application/vnd.ms-excel"
                )
            else:
                st.warning("本次筛选未发现符合条件的 Hits。")
                
        else:
            st.error("分析完成，但没有有效数据（可能全部被剔除或格式错误）。")
            
    except Exception as e:
        st.error(f"分析筛选数据时出错: {e}")
