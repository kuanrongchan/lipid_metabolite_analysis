import pandas as pd
import numpy as np
import io

import seaborn as sns
import matplotlib.pyplot as plt

import streamlit as st
from streamlit_tags import st_tags

def file_upload():
    df_dict, df_names = {}, []
    df_query = st.sidebar.file_uploader("Upload a csv or xlsx file here", type=["csv","xlsx"], accept_multiple_files=False)
    use_demo = st.sidebar.checkbox("Use demo dataset", value=False)
    if df_query is not None:
        head, sep, tail = str(df_query.name).partition(".")
        if tail == 'csv':
            data = st.cache(pd.read_csv)(df_query, index_col=0)
            df_dict[head] = data
            df_names.append(head)

        elif tail == 'xlsx':
            x = st.cache(pd.read_excel)(df_query, index_col=0, sheet_name=None, engine='openpyxl')
            selected_sheet = st.sidebar.multiselect(label="* Select which sheet to read in", options=x.keys())
            for i in selected_sheet:
                df_dict[i] = x[1]
                df_names.append(i)

    else:
        if use_demo:
            data = st.cache(pd.read_csv)("Lipids_VMP1_TMEM41B_KO.csv", index_col=0)
            df_dict["TMEM41B_VMP1_KO"] = data
            df_names.append("TMEM41B_VMP1_KO")
        else:
            st.stop()
    return df_dict, df_names

def log2fc(clean_dict, df_names):
    control_cols = datamanager.multiselect("Select subjects to be used as controls for log2 fold-change calculation", options = clean_dict[df_names[0]].columns.to_list())
    log2FCs = {}
    for k,v in clean_dict.items():
        v = v.select_dtypes(include = np.number)
        log2ct = np.log2(v)
        avg_ct = log2ct.filter(control_cols, axis=1).mean(axis=1)
        avg_ct.name = "average controls"
        logfc = pd.concat([log2ct, avg_ct], axis=1)
        logfc = logfc.subtract(avg_ct, axis=0)
        log2FCs[k] = logfc
    return log2FCs

def data_handler(df_dict, df_names):
    cleaned_dict, metabolites, subjs = {}, [], []
    transpose = datamanager.checkbox("Transpose my dataset", value=False, help="Transpose if your dataset has subjects in the rows and metabolites in columns.")
    for k,v in df_dict.items():
        if transpose:
            numerics = v.select_dtypes(exclude='number')
            df_T = v.T
            df_T = df_T.drop(numerics, axis=0)
            df_T = df_T.astype(float)
            cleaned_dict[k] = df_T
        else:
            df = v.select_dtypes(include='number')
            cleaned_dict[k] = df

    for v in cleaned_dict.values():
        for m in v.index:
            metabolites.append(m)
        for s in v.columns:
            subjs.append(s)

    log2fcs = log2fc(cleaned_dict, df_names)
    filter_metabs = datamanager.multiselect("Select metabolites to filter", options=['all'] + list(set(metabolites)), default='all')
    filter_subjs = datamanager.multiselect("Select subjects to filter", options=['all'] + list(set(subjs)), default='all')

    done = datamanager.checkbox("Finished filtering", value=False)
    if done:
        if 'all' not in filter_metabs:
            for k,v in cleaned_dict.items():
                filtered = v.filter(filter_metabs, axis=0)
                cleaned_dict[k] = filtered
            for k,v in log2fcs.items():
                filtered = v.filter(filter_metabs, axis=0)
                log2fcs[k] = filtered
        else:
            cleaned_dict = cleaned_dict
            log2fcs = log2fcs

        if 'all' not in filter_subjs:
            for k,v in cleaned_dict.items():
                filtered = v.filter(filter_subjs, axis=1)
                cleaned_dict[k] = filtered
            for k,v in log2fcs.items():
                filtered = v.filter(filter_subjs, axis=1)
                log2fcs[k] = filtered
        else:
            cleaned_dict = cleaned_dict
            log2fcs = log2fcs
    else:
        st.stop()
    
    return cleaned_dict, log2fcs

def cluster(df_dict, z_score=None, datatype = 'z-score', width=5, height=10, dendrogram_r = 0.2, dendrogram_c = 0.2, vminmax = (-5,5), cbar_left = 0.98, cbar_bottom = 0.05, cbar_width = 0.15, cbar_height = 0.02):
    if vminmax == (0.0, 0.0):
        vminmax = (None, None)
    for k,df in df_dict.items():
        g = sns.clustermap(df, cmap="vlag",
                        method='average',
                        cbar_pos=(cbar_left, cbar_bottom, cbar_width, cbar_height),
                        center=0, 
                        vmin = vminmax[0], vmax = vminmax[1],
                        z_score=z_score,
                        col_cluster=True,
                        yticklabels=True,
                        figsize=(width, height),
                        dendrogram_ratio=(dendrogram_r, dendrogram_c),
                        linewidths=1, linecolor='white',
                        cbar_kws = {"label":datatype, 'orientation':'horizontal'})

        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=11, rotation=0)
        g.ax_heatmap.set_ylabel("")
        g.fig.suptitle(f"{datatype} metabolites", x=0.5, y=1.02, fontsize=14, fontweight='bold')
        for _, spine in g.ax_heatmap.spines.items():
            spine.set_visible(True)
            spine.set_edgecolor("black")

        st.pyplot(g)
        return g

def create_pdf(fig, fn, graph_module = "pyplot"):
    buf = io.BytesIO()
    if graph_module == 'pyplot':
        fig.savefig(buf, format = 'pdf', bbox_inches = 'tight')
    elif graph_module == 'plotly':
        fig.write_image(file = buf, format = 'pdf')
    st.download_button(
        label = f"Download {fn.replace('_', ' ')} as pdf",
        data = buf,
        file_name = f"{fn}.pdf",
        mime='application/pdf'
        )

######## CONTROL FLOW ###########

st.header("Mass Spectrometry Clustermap")
st.markdown('''

### Instructions to use web tool

1. Upload raw intensity values from mass spectrometry readouts. Data can be uploaded as csv or xlsx files.
2. Ensure that the metabolites are labelled in the first column, and the expression values on the subsequent columns. You may use the web tool to transpose your dataset if the metabolites are in columns.
3. Use the data management expander on the side bar to filter out any metabolites or samples
4. Select the samples that are used as controls. This step is required to plot clustergrams expressed as log2FC values (samples vs controls)
5. Click on the finished filtering checkbox once the settings are finalised.

### Clustermap aesthetics and downloads

- Users can adjust the clustermap dimensions, dendrogram lengths and heights, range of values to be plotted and the colourbar position and size

- You may download the clustergram figures as pdf files by clicking on the download button present at the bottom of the clustergrams.
 ''')

df_dict, df_names = file_upload()
show_df = st.sidebar.checkbox("Show dataset")
datamanager = st.sidebar.expander("Data Management", expanded=True)
if show_df:
    for k, v in df_dict.items():
        st.markdown(f"**{k}**")
        st.info("Click anywhere within the dataframe and use Cmd+F or Ctrl+F to search the dataframe.")
        st.write(v)
df_dict, log2fc_dict = data_handler(df_dict, df_names)

cluster_exp = st.sidebar.expander("Clustermap options", expanded=True)

with cluster_exp:
    width = st.slider("Clustermap width (in inches)", min_value=4, max_value=20, step=1, value = 7, key='width')
    height = st.slider("Clustermap height (in inches)", min_value=4, max_value=40, step=1, value = 14, key='height')
    dendrogram_r = st.slider("Adjust relative row dendrogram length", min_value=0.01, max_value=1.00, step=0.01, value=0.12)
    dendrogram_c = st.slider("Adjust relative column dendrogram height", min_value=0.01, max_value=1.00, step=0.01, value=0.08)
    vminmax = st.slider("Adjust minimum and maximum values of the colour bar", min_value=-10.0, max_value=10.0, step = 0.5, value=(-2.0, 2.0))
    cbar_left = st.slider("Adjust colourbar position from left", min_value=0.0, max_value=1.0, step=0.01, value = 0.78)
    cbar_bottom = st.slider("Adjust colourbar position from bottom", min_value=0.0, max_value=1.0, step=0.01, value = 0.05)
    cbar_width = st.slider("Adjust colourbar width", min_value=0.0, max_value=1.0, step=0.01, value = 0.15)
    cbar_height = st.slider("Adjust colourbar height", min_value=0.0, max_value=1.0, step=0.01, value = 0.02)
    aesthetics_complete = st.checkbox("Finished editing clustermap aesthetics", value=True)


z_cluster, fc_cluster = st.tabs(["Z-score clustermap", 'Log2FC clustermap'])

with z_cluster:
    if aesthetics_complete:
        zfig = cluster(df_dict, z_score=0, datatype='z-score',
        width=width, height=height,
        dendrogram_r = dendrogram_r, dendrogram_c=dendrogram_c,
        vminmax=vminmax,
        cbar_left=cbar_left, cbar_bottom=cbar_bottom, cbar_width=cbar_width, cbar_height=cbar_height
        )
        create_pdf(zfig, "zscore_clustergram", graph_module='pyplot')

with fc_cluster:
    if aesthetics_complete:
        fcfig = cluster(log2fc_dict, datatype='log2FC',
        z_score=None, width=width, height=height,
        dendrogram_r = dendrogram_r, dendrogram_c=dendrogram_c,
        vminmax=vminmax,
        cbar_left=cbar_left, cbar_bottom=cbar_bottom, cbar_width=cbar_width, cbar_height=cbar_height
        )
        create_pdf(fcfig, "log2FC_clustergram", graph_module='pyplot')
