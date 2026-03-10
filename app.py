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
- **支持多板分析**：自动识别板号（每板间隔默认11行）。
- **智能清洗**：自动剔除内参异常值和毒性假阳性。
- **化合物库映射 (New!)**：上传 MCE 等化合物库表格，自动将孔位坐标翻译为具体的药物名称。
""")

# 侧边栏：参数配置
with st.sidebar:
    st.header("⚙️ 参数配置")
    
    PLATE_HEIGHT = st.number_input("每板行数 (Data Rows)", value=8, help="例如 B1-M8 是8行")
    PLATE_STEP = st.number_input("板间距 (Step)", value=11, help="下一板起始行 - 上一板起始行")
    
    st.divider()
    
    THRESHOLD_HIT = st.slider("Hit 判定阈值 (Ratio ≥ 该值)", 1.0, 10.0, 2.0, 0.5)
    LIMIT_HIGH_SIGNAL = st.number_input("高信号异常剔除倍数 (防沉淀/杂质)", value=15.0, help="如果信号高得离谱(比如>15倍)，可能是沉淀，予以剔除")
    
    st.divider()
    st.info("默认配置：\n内参: 第2列 (前5个)\n毒性: 第13列 (前5个)\n药物: 第3-12列")

# 固定列配置
COL_IDX_PLATE_ID = 0   # A列
COL_IDX_CONTROL = 1    # B列
COL_IDX_TOX = 12       # M列
COLS_IDX_DRUG = range(2, 12)
ROWS_CONTROL_RELATIVE = range(0, 5)
ROWS_TOX_RELATIVE = range(0, 5)

# ==========================================
# 1. 核心计算与辅助函数
# ==========================================
def calculate_clean_mean(values):
    """ 计算均值 (剔除 Mean ± 2SD 异常值) """
    data = pd.to_numeric(values, errors='coerce')
    data = data[~np.isnan(data)]
    if len(data) < 3: return np.mean(data) if len(data)>0 else 0
    mean_val = np.mean(data)
    std_val = np.std(data, ddof=1)
    clean = data[(data >= mean_val - 2*std_val) & (data <= mean_val + 2*std_val)]
    return np.mean(clean)

def get_col_letter(col_idx):
    """ Excel 列号转换 (0->A, 1->B...) """
    if 0 <= col_idx <= 25: return chr(65 + col_idx)
    return f"C{col_idx}"

def format_well(w):
    """ 将孔位标准化，例如 A2 -> A02 """
    w = str(w).strip()
    if len(w) >= 2 and w[0].isalpha() and w[1:].isdigit():
        return f"{w[0].upper()}{int(w[1:]):02d}"
    return w

def load_mce_library(file_obj, filename):
    """ 智能读取 MCE 化合物库（自动跳过前几行的介绍文字） """
    try:
        # 先读取前30行寻找表头
        if filename.endswith('.csv'):
            try:
                raw_df = pd.read_csv(file_obj, header=None, nrows=30, encoding='utf-8')
            except UnicodeDecodeError:
                file_obj.seek(0)
                raw_df = pd.read_csv(file_obj, header=None, nrows=30, encoding='gbk')
        else:
            raw_df = pd.read_excel(file_obj, header=None, nrows=30)
            
        header_idx = 0
        for i, row in raw_df.iterrows():
            row_str = " ".join(row.astype(str).fillna("").values).lower()
            # MCE表头通常包含 plate 和 product name 
            if "plate" in row_str and ("well" in row_str or "product name" in row_str):
                header_idx = i
                break
                
        file_obj.seek(0)
        if filename.endswith('.csv'):
            try:
                df = pd.read_csv(file_obj, skiprows=header_idx, encoding='utf-8')
            except:
                file_obj.seek(0)
                df = pd.read_csv(file_obj, skiprows=header_idx, encoding='gbk')
        else:
            df = pd.read_excel(file_obj, skiprows=header_idx)
        return df
    except Exception as e:
        st.error(f"解析化合物库失败: {e}")
        return None

# ==========================================
# 2. 主逻辑
# ==========================================
col_upload1, col_upload2 = st.columns(2)
with col_upload1:
    uploaded_file = st.file_uploader("1️⃣ 必填：请上传 CSV 筛选数据文件", type=["csv"])
with col_upload2:
    uploaded_lib = st.file_uploader("2️⃣ 可选：请上传 MCE 化合物信息库", type=["csv", "xlsx"], help="例如：MCE Library-Detailed Information...xlsx")

if uploaded_file is not None:
    try:
        df = pd.read_csv(uploaded_file, header=None)
        
        all_results = []
        total_tox_drop = 0
        total_high_drop = 0
        global_drug_count = 0
        
        # 进度条
        progress_bar = st.progress(0)
        total_steps = len(range(0, len(df), PLATE_STEP))
        current_step = 0
        
        for start_row in range(0, len(df), PLATE_STEP):
            current_step += 1
            progress_bar.progress(min(current_step / total_steps, 1.0))
            
            end_row = start_row + PLATE_HEIGHT
            if end_row > len(df): break
            
            plate_df = df.iloc[start_row:end_row, :]
            
            # 获取板号
            plate_id = str(plate_df.iloc[0, COL_IDX_PLATE_ID])
            if plate_id == 'nan' or plate_id.strip() == '':
                plate_id = f"Plate_{start_row // PLATE_STEP + 1}"
                
            # 提取纯数字板号 (例如 "Plate_1" -> "1")
            num_match = re.search(r'\d+', plate_id)
            physical_plate = num_match.group(0) if num_match else str(start_row // PLATE_STEP + 1)

            # 计算内参和毒性
            try:
                curr_ctrls = plate_df.iloc[ROWS_CONTROL_RELATIVE, COL_IDX_CONTROL].values
                plate_ctrl_mean = calculate_clean_mean(curr_ctrls)
                
                curr_toxs = plate_df.iloc[ROWS_TOX_RELATIVE, COL_IDX_TOX].values
                plate_tox_threshold = calculate_clean_mean(curr_toxs)
                if plate_ctrl_mean == 0: continue
            except:
                continue

            # 处理药物孔
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
                
                # Excel 相对坐标 (如 C5)
                excel_row_num = start_row + rel_r + 1
                excel_col_char = get_col_letter(abs_c)
                
                # 96孔板物理坐标 (如 A02)
                # rel_r: 0->A, 1->B. abs_c: 2->Col 2 (标准MCE排列: 1是空白/内参, 2-11是药)
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
                st.info("🔄 正在解析并匹配化合物库信息...")
                df_lib = load_mce_library(uploaded_lib, uploaded_lib.name)
                
                if df_lib is not None and "Plate" in df_lib.columns and
