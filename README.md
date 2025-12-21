# NeuroBotanica - Systems Biology Platform for Botanical Therapeutics

**Status:** Phase 1 Core Complete - Dataset & ML Framework  
**Patent Filing:** December 22, 2025  
**Version:** 1.0.0

---

## PROJECT OVERVIEW

NeuroBotanica is an AI-powered systems biology platform that connects traditional botanical knowledge with molecular mechanisms through:

- **Structure-Activity Relationship (SAR) Prediction:** Predict therapeutic effects from molecular structures
- **Dimeric Cannabinoid Discovery:** First-of-its-kind prediction of novel dimeric molecules
- **Multi-Compound Optimization:** Entourage effect quantification and therapeutic formulation
- **FDA Schedule III Documentation:** Regulatory compliance framework for cannabis rescheduling

**Scientific Innovation:**
- 63 cannabinoid compounds with 40+ RDKit molecular descriptors each
- 41 clinical studies from NORML database with condition-based inference
- 120 dimeric cannabinoid predictions with synergy scoring
- Patent-protected prediction methodology (provisional filing Dec 22, 2025)

---

## GETTING STARTED WITH VISUAL STUDIO CODE

### Step 1: Install Required Software

**1. Install Visual Studio Code:**
- Download from: https://code.visualstudio.com/
- Install for your operating system (Windows/Mac/Linux)

**2. Install Python 3.11+:**
- Download from: https://www.python.org/downloads/
- **IMPORTANT:** Check "Add Python to PATH" during installation
- Verify installation: Open terminal and run `python --version`

**3. Install Git (optional but recommended):**
- Download from: https://git-scm.com/downloads
- Used for version control

### Step 2: Open Project in VS Code

**Option A: Download and Extract**
1. Download this entire `neurobotanica_project` folder
2. Open VS Code
3. Click `File > Open Folder`
4. Navigate to the `neurobotanica_project` folder and select it
5. Click `Select Folder` (or `Open` on Mac)

**Option B: Clone from Git (if using version control)**
1. Open VS Code
2. Press `Ctrl+Shift+P` (Windows/Linux) or `Cmd+Shift+P` (Mac)
3. Type "Git: Clone" and press Enter
4. Paste your repository URL
5. Choose a folder location and open when prompted

### Step 3: Install VS Code Extensions (Recommended)

Open Extensions panel (`Ctrl+Shift+X` or `Cmd+Shift+X`) and install:

1. **Python** (by Microsoft) - Essential for Python development
2. **Pylance** (by Microsoft) - Advanced Python language support
3. **Jupyter** (by Microsoft) - If you want to run notebooks
4. **GitLens** (optional) - Enhanced Git capabilities
5. **Better Comments** (optional) - Color-coded comment highlighting

### Step 4: Set Up Python Environment

**1. Open Integrated Terminal in VS Code:**
- Click `Terminal > New Terminal` (or press `Ctrl+` `)
- You should see a terminal at the bottom of VS Code

**2. Create Virtual Environment:**
```bash
# Windows
python -m venv venv
venv\Scripts\activate

# Mac/Linux
python3 -m venv venv
source venv/bin/activate
```

You should see `(venv)` appear at the start of your terminal prompt.

**3. Install Required Packages:**
```bash
pip install --upgrade pip
pip install -r requirements.txt
```

This will install all dependencies:
- `rdkit` - Chemistry calculations
- `scikit-learn` - Machine learning
- `pandas` - Data manipulation
- `numpy` - Numerical computing
- `matplotlib` - Plotting
- `pytest` - Testing framework

**Installation may take 5-10 minutes.** RDKit is a large package.

### Step 5: Verify Installation

**Run the test suite:**
```bash
pytest tests/ -v
```

You should see tests passing (green checkmarks).

**If tests fail:**
- Ensure virtual environment is activated (you see `(venv)` in terminal)
- Ensure all packages installed successfully
- Check Python version: `python --version` (should be 3.11+)

### Step 6: Configure VS Code for Python

**1. Select Python Interpreter:**
- Press `Ctrl+Shift+P` (Windows/Linux) or `Cmd+Shift+P` (Mac)
- Type "Python: Select Interpreter"
- Choose the interpreter from your `venv` folder:
  - Windows: `venv\Scripts\python.exe`
  - Mac/Linux: `venv/bin/python`

**2. Configure Linting (Code Quality):**
- VS Code will prompt to enable linting when you open a `.py` file
- Click "Enable" or select `pylint` as linter
- This shows code quality warnings as you type

**3. Configure Formatting (optional but recommended):**
- Install `black` formatter: `pip install black`
- Open Settings (`Ctrl+,` or `Cmd+,`)
- Search for "Python Formatting Provider"
- Select "black"
- Enable "Format On Save"

### Step 7: Start Coding!

**Open and explore these files:**

1. **`src/analysis/cannabinoid_analyzer.py`** - Load and analyze cannabinoid data
2. **`src/ml_models/sar_predictor.py`** - Train ML models for therapeutic prediction
3. **`src/ml_models/dimeric_predictor.py`** - Predict novel dimeric cannabinoids
4. **`data/training/neurobotanica_complete_dataset_63compounds.json`** - View training data

**Run example analysis:**
```bash
python src/analysis/cannabinoid_analyzer.py
```

---

## PROJECT STRUCTURE

```
neurobotanica_project/
│
├── src/                          # Source code
│   ├── data_extraction/          # NORML citation extraction scripts
│   │   └── norml_parser.py       # Parse clinical studies from NORML
│   │
│   ├── ml_models/                # Machine learning models
│   │   ├── sar_predictor.py      # Structure-activity relationship prediction
│   │   ├── dimeric_predictor.py  # Dimeric cannabinoid screening
│   │   └── entourage_optimizer.py # Multi-compound synergy optimization
│   │
│   └── analysis/                 # Data analysis scripts
│       ├── cannabinoid_analyzer.py    # Analyze cannabinoid properties
│       └── clinical_study_analyzer.py # Analyze clinical study patterns
│
├── data/                         # Data files
│   ├── raw/                      # Raw unprocessed data
│   ├── processed/                # Cleaned and processed data
│   │   └── neurobotanica_descriptors_validated.json
│   │
│   └── training/                 # ML training datasets
│       ├── neurobotanica_complete_dataset_63compounds.json
│       └── neurobotanica_dimeric_predictions.json
│
├── tests/                        # Unit tests and integration tests
│   ├── test_sar_predictor.py
│   ├── test_dimeric_predictor.py
│   └── test_data_loading.py
│
├── docs/                         # Documentation
│   ├── INTEGRATION_GUIDE.md      # Complete technical integration guide
│   ├── API.md                    # API documentation
│   └── PATENT_SUMMARY.md         # Patent claims summary
│
├── requirements.txt              # Python dependencies
├── README.md                     # This file
├── .gitignore                    # Files to exclude from Git
└── pytest.ini                    # Test configuration
```

---

## KEY DATA FILES

### 1. Complete Dataset (63 Cannabinoids)
**File:** `data/training/neurobotanica_complete_dataset_63compounds.json`

**Contents:**
- 12 primary cannabinoids (THC, CBD, CBG, CBC, CBN, THCV, CBDV, etc.)
- 4 active metabolites (11-OH-THC, 7-OH-CBD, THC-COOH, 11-OH-CBD)
- 47 dimeric cannabinoids (predicted novel structures)
- 40+ RDKit molecular descriptors per compound
- Receptor binding profiles (CB1, CB2, 5-HT1A, etc.)
- Therapeutic targets by condition (PTSD, epilepsy, insomnia, Alzheimer's)

**Usage:**
```python
import json

with open('data/training/neurobotanica_complete_dataset_63compounds.json', 'r') as f:
    dataset = json.load(f)

# Access cannabinoid data
thc_data = dataset['compounds']['THC']
print(f"THC Molecular Weight: {thc_data['molecular_weight']}")
print(f"THC LogP: {thc_data['rdkit_descriptors']['logP']}")
print(f"THC CB1 Affinity: {thc_data['receptor_binding']['CB1']['affinity_nM']} nM")
```

### 2. Dimeric Predictions (120 Novel Compounds)
**File:** `data/training/neurobotanica_dimeric_predictions.json`

**Contents:**
- 120 predicted dimeric cannabinoid combinations
- Formation probability scores (0.0-1.0)
- Synergy scores for entourage effects
- Therapeutic potential by condition
- Formation mechanisms (oxidative coupling, Diels-Alder, etc.)

**Usage:**
```python
import json

with open('data/training/neurobotanica_dimeric_predictions.json', 'r') as f:
    dimers = json.load(f)

# Find top dimers for epilepsy
epilepsy_dimers = [
    d for d in dimers['dimeric_cannabinoids'] 
    if d['therapeutic_potential']['Epilepsy'] > 0.8
]
epilepsy_dimers.sort(key=lambda x: x['synergy_score'], reverse=True)

print(f"Top dimer for epilepsy: {epilepsy_dimers[0]['name']}")
print(f"Synergy score: {epilepsy_dimers[0]['synergy_score']}")
```

### 3. Clinical Studies (41 Studies)
**File:** Embedded in `neurobotanica_complete_dataset_63compounds.json` under `clinical_studies` key

**Contents:**
- 41 peer-reviewed clinical studies from NORML database
- Conditions: PTSD (9), Epilepsy (13), Insomnia (10), Alzheimer's (9)
- Cannabinoid profiles (explicit or inferred with confidence weights)
- Study metadata (authors, year, journal, DOI)
- Therapeutic outcomes

---

## COMMON WORKFLOWS

### Workflow 1: Analyze Cannabinoid Properties

```python
# File: src/analysis/cannabinoid_analyzer.py
import json
import pandas as pd

# Load data
with open('data/training/neurobotanica_complete_dataset_63compounds.json', 'r') as f:
    data = json.load(f)

# Convert to DataFrame for analysis
compounds = []
for name, cannabinoid in data['compounds'].items():
    if cannabinoid.get('type') == 'monomeric':  # Exclude dimers for now
        compounds.append({
            'name': name,
            'molecular_weight': cannabinoid['molecular_weight'],
            'logP': cannabinoid['rdkit_descriptors']['logP'],
            'TPSA': cannabinoid['rdkit_descriptors']['TPSA'],
            'psychoactive': cannabinoid.get('psychoactive', False),
            'CB1_affinity': cannabinoid['receptor_binding'].get('CB1', {}).get('affinity_nM')
        })

df = pd.DataFrame(compounds)

# Analysis
print("Cannabinoid Property Statistics:")
print(df.describe())

# Correlation between LogP and psychoactivity
print("\nAverage LogP by psychoactivity:")
print(df.groupby('psychoactive')['logP'].mean())
```

### Workflow 2: Train SAR Prediction Model

```python
# File: src/ml_models/sar_predictor.py
import json
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

# Load training data
with open('data/training/neurobotanica_complete_dataset_63compounds.json', 'r') as f:
    data = json.load(f)

# Prepare features (RDKit descriptors)
X = []
y = []
for study in data['clinical_studies']:
    for cannabinoid_name in study['cannabinoid_profile']:
        cannabinoid = data['compounds'][cannabinoid_name]
        
        # Extract RDKit descriptors as features
        descriptors = cannabinoid['rdkit_descriptors']
        features = [
            descriptors['logP'],
            descriptors['TPSA'],
            descriptors['molecular_weight'],
            descriptors['h_bond_donors'],
            descriptors['h_bond_acceptors'],
            descriptors['rotatable_bonds'],
            descriptors['num_aromatic_rings']
        ]
        
        X.append(features)
        
        # Label: 1 if effective for condition, 0 otherwise
        y.append(1 if study['outcome'] == 'positive' else 0)

X = np.array(X)
y = np.array(y)

# Train model
model = RandomForestClassifier(n_estimators=100, random_state=42)
scores = cross_val_score(model, X, y, cv=5)

print(f"Cross-validation accuracy: {scores.mean():.2f} (+/- {scores.std():.2f})")

# Train on full dataset
model.fit(X, y)
```

### Workflow 3: Screen Dimeric Compounds

```python
# File: src/ml_models/dimeric_predictor.py
import json

# Load dimeric predictions
with open('data/training/neurobotanica_dimeric_predictions.json', 'r') as f:
    dimers = json.load(f)

# Screen for high-potential dimers
condition = "Epilepsy"
min_synergy = 0.7
min_formation = 0.5

candidates = [
    dimer for dimer in dimers['dimeric_cannabinoids']
    if (dimer['synergy_score'] >= min_synergy and 
        dimer['formation_probability'] >= min_formation and
        dimer['therapeutic_potential'].get(condition, 0) >= 0.8)
]

# Sort by composite score
candidates.sort(
    key=lambda x: (
        x['synergy_score'] * 0.4 + 
        x['formation_probability'] * 0.3 +
        x['therapeutic_potential'].get(condition, 0) * 0.3
    ),
    reverse=True
)

print(f"Top 5 dimeric candidates for {condition}:")
for i, dimer in enumerate(candidates[:5], 1):
    print(f"{i}. {dimer['name']}")
    print(f"   Synergy: {dimer['synergy_score']:.2f}")
    print(f"   Formation Probability: {dimer['formation_probability']:.2f}")
    print(f"   Therapeutic Potential: {dimer['therapeutic_potential'][condition]:.2f}")
    print(f"   Formation Mechanism: {dimer['formation_mechanism']}")
    print()
```

---

## RUNNING TESTS

**Run all tests:**
```bash
pytest tests/ -v
```

**Run specific test file:**
```bash
pytest tests/test_sar_predictor.py -v
```

**Run with coverage report:**
```bash
pytest tests/ --cov=src --cov-report=html
```

Then open `htmlcov/index.html` in your browser to see detailed coverage.

---

## DEVELOPMENT WORKFLOW

### 1. Daily Development
1. Activate virtual environment: `source venv/bin/activate` (or `venv\Scripts\activate` on Windows)
2. Pull latest changes (if using Git): `git pull`
3. Run tests to ensure everything works: `pytest tests/ -v`
4. Make your changes to code
5. Run tests again to verify changes: `pytest tests/ -v`
6. Commit changes (if using Git): `git add .` then `git commit -m "Description of changes"`

### 2. Adding New Features
1. Create new file in appropriate `src/` subdirectory
2. Write function with docstring and type hints
3. Create corresponding test in `tests/`
4. Run test: `pytest tests/test_your_new_feature.py -v`
5. Iterate until test passes
6. Update documentation in `docs/`

### 3. Data Updates
1. Place new raw data in `data/raw/`
2. Process data using scripts in `src/data_extraction/`
3. Save processed data to `data/processed/`
4. Update training dataset in `data/training/`
5. Retrain models using scripts in `src/ml_models/`
6. Validate model performance with new data

---

## TROUBLESHOOTING

### "Module not found" Error
**Problem:** Python can't find your modules
**Solution:** 
1. Ensure virtual environment is activated (`(venv)` in terminal prompt)
2. Ensure you're in the project root directory (`neurobotanica_project/`)
3. Install requirements again: `pip install -r requirements.txt`

### RDKit Installation Fails
**Problem:** RDKit is difficult to install on some systems
**Solution:**
- **Recommended:** Use Anaconda instead of pip
  ```bash
  conda create -n neurobotanica python=3.11
  conda activate neurobotanica
  conda install -c conda-forge rdkit
  pip install -r requirements.txt  # Install remaining packages
  ```

### VS Code Not Recognizing Virtual Environment
**Problem:** Autocomplete not working, imports show errors
**Solution:**
1. Press `Ctrl+Shift+P` (or `Cmd+Shift+P`)
2. Type "Python: Select Interpreter"
3. Choose the interpreter from your `venv` folder
4. Reload VS Code window: `Ctrl+Shift+P` → "Developer: Reload Window"

### Tests Fail with "No such file or directory"
**Problem:** Tests can't find data files
**Solution:**
1. Ensure you're running tests from project root: `cd neurobotanica_project/`
2. Verify data files exist: `ls data/training/`
3. Check file paths in test files match actual locations

---

## NEXT STEPS

### For Development:
1. ✅ **Phase 1 Complete:** Dataset ready, models trainable
2. ⏳ **Phase 2 (Months 1-3):** Expand dataset to 200 studies
3. 📅 **Phase 3 (Months 4-6):** Integrate terpene profiles, expand to 500 studies
4. 📅 **Phase 4 (Months 7-12):** Production deployment, continuous learning pipeline

### For Learning:
1. Read `docs/INTEGRATION_GUIDE.md` for complete technical details
2. Explore data files in `data/` to understand structure
3. Run example scripts in `src/analysis/`
4. Modify scripts to answer your own research questions
5. Read scientific papers on cannabinoid pharmacology

### For Collaboration:
1. Use Git for version control (see Git basics below)
2. Create feature branches for new work
3. Write tests for all new code
4. Document your code with docstrings
5. Submit pull requests for review

---

## GIT BASICS (OPTIONAL)

**Initialize repository:**
```bash
git init
git add .
git commit -m "Initial commit - Phase 1 complete"
```

**Create feature branch:**
```bash
git checkout -b feature/terpene-integration
# Make changes
git add .
git commit -m "Add terpene molecular descriptors"
git checkout main
git merge feature/terpene-integration
```

**Connect to remote repository:**
```bash
git remote add origin https://github.com/yourusername/neurobotanica.git
git push -u origin main
```

---

## RESOURCES

### Documentation:
- Full integration guide: `docs/INTEGRATION_GUIDE.md`
- API documentation: `docs/API.md`
- Patent summary: `docs/PATENT_SUMMARY.md`

### Python Libraries:
- RDKit documentation: https://www.rdkit.org/docs/
- scikit-learn tutorials: https://scikit-learn.org/stable/tutorial/
- pandas guide: https://pandas.pydata.org/docs/user_guide/

### VS Code:
- Python in VS Code: https://code.visualstudio.com/docs/python/python-tutorial
- Testing in VS Code: https://code.visualstudio.com/docs/python/testing

### Domain Knowledge:
- NORML clinical studies: https://norml.org/marijuana/library/recent-research-on-medical-marijuana/
- Cannabis pharmacology: PubMed searches for "cannabinoid pharmacology"
- Systems biology: Nature Systems Biology journal

---

## SUPPORT

**For technical issues:**
- Check troubleshooting section above
- Review error messages carefully
- Search Stack Overflow for specific errors
- Check Python/RDKit documentation

**For scientific questions:**
- Review clinical studies in dataset
- Read peer-reviewed literature on cannabinoids
- Consult systems biology textbooks

---

## LICENSE & PATENT

**Software License:** Proprietary (Cloak and Quill Research 501(c)(3))  
**Patent Status:** Provisional application filed December 22, 2025  
**Trade Secrets:** Analytical processing modules (not included in this release)

**For commercial licensing inquiries, contact:** [Your Contact Information]

---

**Last Updated:** December 16, 2025  
**Version:** 1.0.0  
**Phase 1 Status:** ✅ Complete
