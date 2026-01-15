Guidance for the Reproducibility Reports

**Data Files**:

**Synthetic Data**: 
The synthetic datasets are simulated through the code within **FS_CI_Final.R** / **FS_CI_SL.R**, with various settings.

**Final_OUD.csv & Final_NO_OUD.csv**: 
These two datasets are derived from the National Survey on Drug Use and Health (NSDUH).
- Data Source: https://www.samhsa.gov/data/data-we-collect/nsduh-national-survey-drug-use-and-health
- Access Instructions: Navigate to "Download Data Files" section, select the 2015-2019 year, and download the dataset. The raw NSDUH data should be processed using **nsduh_update.R** to generate Final_OUD.csv and Final_NO_OUD.csv.

**diab.csv**: 
CDC Diabetes Health Indicators dataset from UCI Machine Learning Repository.
- Direct Link: https://archive.ics.uci.edu/dataset/891/cdc+diabetes+health+indicators
- Access Instructions: Download the dataset from the "Download" section on the dataset page. Rename the file to **diab.csv**.

**obs.csv**: 
Estimation of Obesity Levels Based on Eating Habits and Physical Condition dataset from UCI Machine Learning Repository.
- Direct Link: https://archive.ics.uci.edu/dataset/544/estimation+of+obesity+levels+based+on+eating+habits+and+physical+condition
- Access Instructions: Download the dataset from the "Download" section on the dataset page. Rename the file to **obs.csv**.

**Code Files**:

**FS_CI_Final.R:** This file compares the existing state-of-the-art feature selection frameworks with the proposed three-stage framework using simulated datasets.

**FS_CI_SL.R:** This file compares the performance before and after feature selection with the proposed three-stage framework.

**sensitivity.R:** Sensitivity analysis for hyperparameters (gamma_1 and gamma_2) for the three-stage feature selection framework.

**nsduh_update.R:** Implement three-stage feature selection framework to the Opioid-Use dataset (Final_NO_OUD.csv/Final_OUD.csv). This script also contains data preprocessing steps.

**obs.R:** Implement three-stage feature selection to the EOL-EHPC dataset (obs.csv).

**dia.R:** Implement three-stage feature selection to the diabetes dataset (diab.csv).

**Other Notes**:

**Required R Package - lqa**:
The lqa package is required for reproducing the OAL method (Shortreed. et al. 2017).
- Package Source: https://cran.r-project.org/src/contrib/Archive/lqa/
- Installation: Download and install manually using install.packages("path/to/lqa_x.x-x.tar.gz", repos = NULL, type = "source")

**Contact Information**:
For questions regarding data access or code reproduction:
- Tianyu Yang: yang.tianyu@northeastern.edu
- Md. Noor-E-Alam: mnalam@neu.edu
