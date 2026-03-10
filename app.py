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
- **支持多板分析**：自动识别板号。
- **智能清洗**：自动剔除内参异常值和毒性假阳性。
- **化合物库映射**：支持 MCE 的各种刁钻格式 (CSV/TSV/Excel)，自动匹配药物名！
""")

# 侧边栏：参数配置
with st.sidebar:
    st.header("⚙️ 参数配置")
    PLATE_HEIGHT = st.number_input("每板行数 (Data Rows)", value=8, help="例如 B1-M8 是8行")
    PLATE_STEP = st.number_input("板间距 (Step)", value=11, help="下一板起始行 - 上一板起始行")
    st.divider()
    THRESHOLD_HIT = st.slider("Hit 判定阈值 (Ratio ≥ 该值)", 1.0, 10.0, 2.0, 0.5)
    LIMIT_HIGH_SIGNAL = st.number_input("高信号异常剔除倍数 (防沉淀/杂质)", value=15.0)
    st.divider()
    st.info("默认配置：\n内参: 第2列\n毒性: 第13列\n药物: 第3-12列")

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
    """ 彻底清洗孔位：删掉空格、Tab、换行符 """
    w = str(w)
    w = re.sub(r'[\s\t\n\r\ufeff]+', '', w)
    if len(w) >= 2 and w[0].isalpha() and w[1:].isdigit():
        return f"{w[0].upper()}{int(w[1:]):02d}"
    return w

def load_mce_library(file_obj, filename):
    """ 超级健壮的化合物库读取器（自动探测逗号或Tab） """
    try:
        if filename.endswith('.csv'):
            # 1. 强行读取文件内容为字符串
            file_obj.seek(0)
            raw_bytes = file_obj.read()
            try:
                text = raw_bytes.decode('utf-8-sig') # 自动去除UTF-8 BOM
            except UnicodeDecodeError:
                text = raw_bytes.decode('gbk')
            
            lines = text.splitlines()
            header_idx = 0
            sep = ',' # 默认用逗号
            
            # 2. 扫描前30行找表头，并智能判断分隔符
            for i, line in enumerate(lines[:30]):
                line_lower = line.lower()
                if 'plate' in line_lower and ('well' in line_lower or 'product' in line_lower or 'catalog' in line_lower):
                    header_idx = i
                    # 如果这行里有 \t，说明它其实是 TSV
                    if '\t' in line:
                        sep = '\t'
                    break
            
            # 3. 使用找到的分隔符去解析字符串
            df = pd.read_csv(io.StringIO(text), skiprows=header_idx, sep=sep, engine='python')
            
        else:
            # Excel 处理逻辑
            df = pd.read_excel(file_obj)
            cols_str = " ".join(df.columns.astype(str)).lower()
            if "plate" not in cols_str:
                header_idx = -1
                for i, row in df.head(30).iterrows():
                    row_str = " ".join(row.astype(str).fillna("")).lower()
                    if "plate" in row_str and ("well" in row_str or "product" in row_str):
                        header_idx = i
                        break
                if header_idx != -1:
                    df.columns = df.iloc[header_idx]
                    df = df.iloc[header_idx+1:].reset_index(drop=True)

        # 4. 标准化列名，彻底清理那些看不见的牛鬼蛇神符号
        df.columns = df.columns.astype(str).str.strip().str.replace('\ufeff', '')
        
        col_map = {}
        for c in df.columns:
            cl = c.lower().strip()
            if cl in ['plate', 'plate number', 'platenum']: col_map[c] = 'Plate'
            elif cl in ['well', 'position', 'well position']: col_map[c] = 'Well'
            elif 'product name' in cl or 'compound name' in cl or cl == 'name': col_map[c] = 'Product Name'
            elif 'catalog number' in cl or 'catalog' in cl or 'item' in cl: col_map[c] = 'Catalog Number'
            elif 'target' in cl: col_map[c] = 'Target'
            elif 'pathway' in cl: col_map[c] = 'Pathway'
        
        df.rename(columns=col_map, inplace=True)
        return df
    except Exception as e:
        st.error(f"解析化合物库失败: {e}")
        return None

# ==========================================
# 2. 主逻辑
# ==========================================
col_upload1, col_upload2 = st.columns(2)
with col_upload1:
    uploaded_file = st.file_uploader("1️⃣ 必填：请上传 CSV 筛选数据", type=["csv"])
with col_upload2:
    uploaded_lib = st.file_uploader("2️⃣ 可选：请上传 MCE 化合物库 (CSV/Excel)", type=["csv", "xlsx"])

if uploaded_file is not None:
    try:
        df = pd.read_csv(uploaded_file, header=None)
        
        all_results = []
        total_tox_drop = 0
        total_high_drop = 0
        global_drug_count = 0
        
        progress_bar = st.progress(0)
        total_steps = len(range(0, len(df), PLATE_STEP))
        current_step = 0
        
        for start_row in range(0, len(df), PLATE_STEP):
            current_step += 1
            progress_bar.progress(min(current_step / total_steps, 1.0))
            
            end_row = start_row + PLATE_HEIGHT
            if end_row > len(df): break
            
            plate_df = df.iloc[start_row:end_row, :]
            
            plate_id = str(plate_df.iloc[0, COL_IDX_PLATE_ID])
            if plate_id == 'nan' or plate_id.strip() == '':
                plate_id = f"Plate_{start_row // PLATE_STEP + 1}"
                
            num_match = re.search(r'\d+', plate_id)
            physical_plate = num_match.group(0) if num_match else str(start_row // PLATE_STEP + 1)

            try:
                curr_ctrls = plate_df.iloc[ROWS_CONTROL_RELATIVE, COL_IDX_CONTROL].values
                plate_ctrl_mean = calculate_clean_mean(curr_ctrls)
                
                curr_toxs = plate_df.iloc[ROWS_TOX_RELATIVE, COL_IDX_TOX].values
                plate_tox_threshold = calculate_clean_mean(curr_toxs)
                if plate_ctrl_mean == 0: continue
            except:
                continue

            drug_block = plate_df.iloc[:, COLS_IDX_DRUG]
            
            for (rel_r, abs_c), val in drug_block.stack().items():
                global_drug_count += 1
                try: val = float(val)
                except: continue
                if np.isnan(val): continue

                if val < plate_tox_threshold:
                    total_tox_drop += 1
                    continue
                if val > (plate_ctrl_mean * LIMIT_HIGH_SIGNAL):
                    total_high_drop += 1
                    continue
                
                ratio = val / plate_ctrl_mean
                is_hit = "Yes" if ratio >= THRESHOLD_HIT else "No"
                
                excel_row_num = start_row + rel_r + 1
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
                    "Raw_RFU": val,
                    "Is_Hit": is_hit
                })

        # ==========================
        # 合并化合物信息库
        # ==========================
        if all_results:
            df_res = pd.DataFrame(all_results)
            
            if uploaded_lib is not None:
                st.info("🔄 正在智能解析化合物库...")
                df_lib = load_mce_library(uploaded_lib, uploaded_lib.name)
                
                if df_lib is not None:
                    if "Plate" not in df_lib.columns or "Well" not in df_lib.columns:
                        st.warning("⚠️ 格式解析异常，找不到标准列名。已提取的列名如下：")
                        st.write(df_lib.columns.tolist())
                    else:
                        # 确保 Plate 和 Well 格式统一
                        df_lib["Physical_Plate"] = df_lib["Plate"].astype(str).str.extract(r'(\d+)')[0]
                        df_lib["Physical_Well"] = df_lib["Well"].apply(format_well)
                        
                        df_res = pd.merge(df_res, df_lib, on=["Physical_Plate", "Physical_Well"], how="left")
                        st.success("✅ 化合物名称匹配成功！")

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
                # 重新组织显示的列，把化合物名称放在最前面
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
                    label="📥 下载详细 Hits 列表 (Excel)",
                    data=buffer.getvalue(),
                    file_name="Enhancer_Hits_with_Info.xlsx",
                    mime="application/vnd.ms-excel"
                )
            else:
                st.warning("本次筛选未发现符合条件的 Hits。")
                
        else:
            st.error("分析完成，但没有有效数据（可能全部被剔除或格式错误）。")
            
    except Exception as e:
        st.error(f"文件读取或分析出错: {e}")
