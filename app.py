import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import io

# ==========================================
# 0. 网页配置与参数设置
# ==========================================
st.set_page_config(page_title="FeliEase 药物筛选分析平台 (激活剂模式)", layout="wide")

st.title("📈 FeliEase 高通量药物筛选分析平台 - 激活剂/升高信号筛选")
st.markdown("""
此工具用于自动化分析药物筛选数据，**寻找引起信号升高的化合物**。
- **支持多板分析**：自动识别板号（每板间隔默认11行）。
- **智能清洗**：自动剔除内参异常值和毒性假阳性。
- **自动归一化**：基于每板内参计算 Ratio，自动标记高信号 Hit。
""")

# 侧边栏：参数配置
with st.sidebar:
    st.header("⚙️ 参数配置")
    
    PLATE_HEIGHT = st.number_input("每板行数 (Data Rows)", value=8, help="例如 B1-M8 是8行")
    PLATE_STEP = st.number_input("板间距 (Step)", value=11, help="下一板起始行 - 上一板起始行")
    
    st.divider()
    
    # 【修改点 1】滑块范围调整，默认寻找 > 2.0 的药
    THRESHOLD_HIT = st.slider("Hit 判定阈值 (Ratio ≥ 该值)", 1.0, 10.0, 2.0, 0.5)
    
    # 【修改点 2】把异常值剔除上限调高，防止误删强效高信号药
    LIMIT_HIGH_SIGNAL = st.number_input("高信号异常剔除倍数 (防沉淀/杂质)", value=15.0, help="如果信号高得离谱(比如>15倍)，可能是化合物自身发光或沉淀，予以剔除")
    
    st.divider()
    st.info("默认配置：\n内参: 第2列 (前5个)\n毒性: 第13列 (前5个)\n药物: 第3-12列")

# 固定列配置 (如需修改请直接改代码)
COL_IDX_PLATE_ID = 0   # A列
COL_IDX_CONTROL = 1    # B列
COL_IDX_TOX = 12       # M列
COLS_IDX_DRUG = range(2, 12)
ROWS_CONTROL_RELATIVE = range(0, 5)
ROWS_TOX_RELATIVE = range(0, 5)

# ==========================================
# 1. 核心计算函数
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
    if 0 <= col_idx <= 25: return chr(65 + col_idx)
    return f"C{col_idx}"

# ==========================================
# 2. 主逻辑
# ==========================================
uploaded_file = st.file_uploader("📂 请上传 CSV 数据文件", type=["csv"])

if uploaded_file is not None:
    try:
        # 读取数据
        df = pd.read_csv(uploaded_file, header=None)
        st.success(f"✅ 文件读取成功！共 {len(df)} 行数据。")
        
        all_results = []
        total_tox_drop = 0
        total_high_drop = 0
        global_drug_count = 0
        
        # 进度条
        progress_bar = st.progress(0)
        
        # 循环处理每一板
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

                # 判定逻辑
                if val < plate_tox_threshold:
                    total_tox_drop += 1
                    continue
                
                if val > (plate_ctrl_mean * LIMIT_HIGH_SIGNAL):
                    total_high_drop += 1
                    continue
                
                ratio = val / plate_ctrl_mean
                
                # 【修改点 3】Hit 判定逻辑翻转：大于等于阈值才是 Hit
                is_hit = "Yes" if ratio >= THRESHOLD_HIT else "No"
                
                excel_row_num = start_row + rel_r + 1
                excel_col_char = get_col_letter(abs_c)
                
                all_results.append({
                    "Global_ID": global_drug_count,
                    "Plate_ID": plate_id,
                    "Coordinate": f"{excel_col_char}{excel_row_num}",
                    "Ratio": ratio,
                    "Raw_RFU": val,
                    "Is_Hit": is_hit,
                    "Excel_Row": excel_row_num,
                    "Excel_Col": excel_col_char
                })

        # ==========================
        # 结果展示
        # ==========================
        if all_results:
            df_res = pd.DataFrame(all_results)
            
            # 1. 统计信息
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("总扫描孔数", global_drug_count)
            col2.metric("毒性/假阳性剔除", total_tox_drop)
            col3.metric("异常极高信号剔除", total_high_drop)
            col4.metric("有效数据量", len(df_res))
            
            st.divider()
            
            # 2. 绘图
            st.subheader("📊 筛选结果可视化")
            fig, ax = plt.subplots(figsize=(12, 6))
            
            # 普通点 (绿色)
            ax.scatter(df_res["Global_ID"], df_res["Ratio"], c='#5cb85c', s=15, alpha=0.6, label='Candidates')
            
            # Hits 点 (红色)
            hits_df = df_res[df_res["Is_Hit"] == "Yes"]
            if not hits_df.empty:
                ax.scatter(hits_df["Global_ID"], hits_df["Ratio"], c='red', s=40, label=f'Hits (≥ {THRESHOLD_HIT})')
                
            ax.axhline(1.0, color='black', linestyle='--', label='Plate Control (1.0)')
            
            # 画一条表示 Hit 阈值的参考线 (蓝色虚线)
            ax.axhline(THRESHOLD_HIT, color='blue', linestyle=':', label=f'Hit Threshold ({THRESHOLD_HIT})')
            
            ax.set_xlabel('Global Drug Sequence')
            ax.set_ylabel('Normalized Ratio')
            
            # 调整Y轴，确保能看全所有高信号点
            y_max = df_res["Ratio"].max()
            ax.set_ylim(bottom=-0.1, top=max(y_max * 1.1, THRESHOLD_HIT * 1.5))
            ax.legend()
            ax.grid(True, linestyle=':', alpha=0.4)
            
            st.pyplot(fig)
            
            # 3. Hits 下载
            st.divider()
            st.subheader("🎉 强效升高信号药 (Hits) 列表")
            
            if not hits_df.empty:
                st.write(f"发现 **{len(hits_df)}** 个 Hits (Ratio ≥ {THRESHOLD_HIT})")
                
                # 【修改点 4】按 Ratio 降序排列 (ascending=False)，最高的排在最前面
                hits_display = hits_df.sort_values(by="Ratio", ascending=False)[["Plate_ID", "Coordinate", "Ratio", "Raw_RFU", "Global_ID"]]
                st.dataframe(hits_display.head(10))
                
                # 下载按钮 (Excel)
                buffer = io.BytesIO()
                with pd.ExcelWriter(buffer, engine='xlsxwriter') as writer:
                    hits_display.to_excel(writer, index=False, sheet_name='Hits')
                
                st.download_button(
                    label="📥 下载 Hits 列表 (Excel)",
                    data=buffer.getvalue(),
                    file_name="Enhancer_Hits.xlsx",
                    mime="application/vnd.ms-excel"
                )
            else:
                st.warning("本次筛选未发现符合条件的 Hits。")
                
        else:
            st.error("分析完成，但没有有效数据（可能全部被剔除或格式错误）。")
            
    except Exception as e:
        st.error(f"文件读取或分析出错: {e}")
