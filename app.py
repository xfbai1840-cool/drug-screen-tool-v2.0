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
- **自动定位板位**：无需设置板间距，程序自动识别每板的起始位置。
- **智能清洗**：自动剔除内参异常值和毒性假阳性。
- **化合物映射**：请上传 **MCE 的 Excel (.xlsx) 信息表**，程序会自动剔除干扰行并匹配药物名称！
""")

# 侧边栏：参数配置
with st.sidebar:
    st.header("⚙️ 参数配置")
    st.success("🤖 自动找板引擎已激活\n无需手动设置板间距！")
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
    """ 暴力清洗孔位：删掉空格、Tab、换行符 """
    w = str(w)
    w = re.sub(r'[\s\t\n\r\ufeff]+', '', w)
    if len(w) >= 2 and w[0].isalpha() and w[1:].isdigit():
        return f"{w[0].upper()}{int(w[1:]):02d}"
    return w

def load_mce_library(file_obj, filename):
    """ 专门针对 MCE Excel 文件的超级解析器 """
    try:
        # 1. 读取 Excel
        df = pd.read_excel(file_obj)
        
        # 2. 寻找真实表头
        header_idx = -1
        for i in range(min(30, len(df))):
            row_str = " ".join(df.iloc[i].fillna("").astype(str)).lower()
            if "plate" in row_str and "well" in row_str:
                header_idx = i
                break
                
        if header_idx != -1:
            df.columns = df.iloc[header_idx].astype(str).str.strip().str.replace('\ufeff', '')
            df = df.iloc[header_idx+1:].reset_index(drop=True)
        else:
            st.warning("⚠️ 无法在 Excel 中找到 Plate 和 Well 表头。")
            return None
            
        # 3. 规范化列名
        col_map = {}
        for c in df.columns:
            cl = c.lower().strip()
            if cl in ['plate', 'rack']: col_map[c] = 'Plate'
            elif cl in ['well', 'position']: col_map[c] = 'Well'
            elif 'product name' in cl or 'compound name' in cl or cl == 'name': col_map[c] = 'Product Name'
            elif 'catalog' in cl or 'item' in cl: col_map[c] = 'Catalog Number'
            elif 'target' in cl: col_map[c] = 'Target'
        df.rename(columns=col_map, inplace=True)
        
        # 4. 【核心】干掉带有 10mM, μL, DMSO 的垃圾注释行
        # 标准孔位是 A01, B02 这种格式，正则匹配 ^[A-Za-z]\d{1,2}$
        if 'Well' in df.columns:
            valid_well_mask = df['Well'].fillna("").astype(str).str.strip().str.match(r'^[a-zA-Z]\d{1,2}$')
            df = df[valid_well_mask].copy()
            
        # 5. 【核心】把 MCE 的乱码板号映射为 1, 2, 3...
        if 'Plate' in df.columns:
            unique_plates = df['Plate'].dropna().unique()
            # 过滤空字符串
            unique_plates = [p for p in unique_plates if p.strip() != '']
            # 建立映射字典：第一个出现的板=1，第二个=2...
            plate_mapping = {p: str(i+1) for i, p in enumerate(unique_plates)}
            df['Physical_Plate'] = df['Plate'].map(plate_mapping)
            
        if 'Well' in df.columns:
            df['Physical_Well'] = df['Well'].apply(format_well)
            
        return df
    except Exception as e:
        st.error(f"解析 Excel 化合物库失败: {e}")
        return None

# ==========================================
# 2. 主逻辑
# ==========================================
col_upload1, col_upload2 = st.columns(2)
with col_upload1:
    uploaded_file = st.file_uploader("1️⃣ 必填：请上传 CSV 筛药结果", type=["csv"])
with col_upload2:
    uploaded_lib = st.file_uploader("2️⃣ 可选：请上传 MCE 库文件 (强烈建议传 .xlsx)", type=["csv", "xlsx"])

if uploaded_file is not None:
    try:
        # 读取筛药结果数据
        df = pd.read_csv(uploaded_file, header=None)
        
        all_results = []
        total_tox_drop = 0
        total_high_drop = 0
        global_drug_count = 0
        plate_count = 1  # 记录物理板数
        
        progress_bar = st.progress(0)
        
        # ---------------------------------------------
        # 【超级找板引擎】自动无视空行，扫描 A 列识别板
        # ---------------------------------------------
        i = 0
        while i < len(df):
            val = str(df.iloc[i, COL_IDX_PLATE_ID]).strip()
            
            # 如果 A 列有内容，且不是 nan，说明找到了一块新板
            if val and val.lower() != 'nan':
                # 检查剩下的行数够不够一整板 (默认 8 行)
                if i + int(PLATE_HEIGHT) <= len(df):
                    plate_df = df.iloc[i : i + int(PLATE_HEIGHT), :].reset_index(drop=True)
                    plate_id = val
                    physical_plate = str(plate_count) # 强制赋予物理序号 1, 2, 3...
                    
                    # --- 计算内参和毒性 ---
                    try:
                        curr_ctrls = plate_df.iloc[ROWS_CONTROL_RELATIVE, COL_IDX_CONTROL].values
                        plate_ctrl_mean = calculate_clean_mean(curr_ctrls)
                        
                        curr_toxs = plate_df.iloc[ROWS_TOX_RELATIVE, COL_IDX_TOX].values
                        plate_tox_threshold = calculate_clean_mean(curr_toxs)
                    except:
                        plate_ctrl_mean = 0
                        
                    if plate_ctrl_mean > 0:
                        # --- 处理药物孔 ---
                        drug_block = plate_df.iloc[:, COLS_IDX_DRUG]
                        for (rel_r, abs_c), cell_val in drug_block.stack().items():
                            global_drug_count += 1
                            try: 
                                cell_val = float(cell_val)
                            except: 
                                continue
                                
                            if np.isnan(cell_val): continue

                            # 判定假阳性和超高信号
                            if cell_val < plate_tox_threshold:
                                total_tox_drop += 1
                                continue
                            if cell_val > (plate_ctrl_mean * LIMIT_HIGH_SIGNAL):
                                total_high_drop += 1
                                continue
                            
                            ratio = cell_val / plate_ctrl_mean
                            is_hit = "Yes" if ratio >= THRESHOLD_HIT else "No"
                            
                            # 计算坐标
                            excel_row_num = i + rel_r + 1
                            excel_col_char = get_col_letter(abs_c)
                            
                            plate_row_char = chr(65 + rel_r)
                            physical_well = f"{plate_row_char}{abs_c:02d}"
                            
                            all_results.append({
                                "Global_ID": global_drug_count,
                                "Plate_ID": plate_id,            # 原始名称，比如 1#-F
                                "Physical_Plate": physical_plate, # 物理板号，比如 1
                                "Physical_Well": physical_well,  # 物理孔位，比如 A02
                                "Coordinate": f"{excel_col_char}{excel_row_num}",
                                "Ratio": ratio,
                                "Raw_RFU": cell_val,
                                "Is_Hit": is_hit
                            })
                    
                    plate_count += 1
                    i += int(PLATE_HEIGHT) # 跳过这 8 行
                    progress_bar.progress(min(i / len(df), 1.0))
                else:
                    break # 剩余行数不够一板，退出
            else:
                i += 1 # 如果是空行，继续往下找

        progress_bar.progress(1.0)

        # ==========================
        # 合并化合物信息库
        # ==========================
        if all_results:
            df_res = pd.DataFrame(all_results)
            
            if uploaded_lib is not None:
                st.info("🔄 正在智能解析 MCE Excel 化合物库...")
                df_lib = load_mce_library(uploaded_lib, uploaded_lib.name)
                
                if df_lib is not None:
                    if "Physical_Plate" not in df_lib.columns or "Physical_Well" not in df_lib.columns:
                        st.warning("⚠️ Excel 格式不匹配，未能生成对照列。")
                    else:
                        # 强制转换为字符串后进行完美合并
                        df_res['Physical_Plate'] = df_res['Physical_Plate'].astype(str)
                        df_res['Physical_Well'] = df_res['Physical_Well'].astype(str)
                        df_lib['Physical_Plate'] = df_lib['Physical_Plate'].astype(str)
                        df_lib['Physical_Well'] = df_lib['Physical_Well'].astype(str)
                        
                        df_res = pd.merge(df_res, df_lib, on=["Physical_Plate", "Physical_Well"], how="left")
                        st.success("✅ 化合物名称已成功接入！")

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
                # 定制展示列：把 MCE 的关键信息排在最前
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
