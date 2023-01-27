# Lipid metabolite analysis webtool

This webtool is built to facilitate researchers visualise their metabolite data quickly. To use this webtool, users will have to:

1. Upload mass spectrometry csv or xlsx files (without log-transformation)
2. Tick the checkbox to transpose your dataset if your subjects are in rows and metabolites are in columns
3. Use the data manager expander to filter out any metabolites or subjects
4. Select the subjects to be used as controls for log2 fold-change calculation (subject vs averaged control)
5. Tick the finished filtering checkbox once filtering is complete.

# Running this webtool locally

## Technical Requirements
Please install the following:

- Install Python 3.7 or later (!=3.9.7) at https://www.python.org/downloads/

## Installing Streamlit locally

```
pip install streamlit
```

## Download this repository

- Click the green 'Code' button on the top right of the repository
- Download ZIP
- Unzip this folder

## Set-up on command line

```
cd path/to/repository
pip install -r requirements.txt
streamlit run app.py
```
