# QUICK START GUIDE - Visual Studio Code Setup

## 🚀 Get Started in 5 Minutes

### Step 1: Download Your Project
1. Download the entire `neurobotanica_project` folder from the outputs
2. Save it somewhere easy to find (e.g., `Documents/neurobotanica_project`)

### Step 2: Install Software (One-Time Setup)

**A. Install Visual Studio Code**
- Go to: https://code.visualstudio.com/
- Click "Download" for your operating system
- Install and open VS Code

**B. Install Python 3.11+**
- Go to: https://www.python.org/downloads/
- Download Python 3.11 or newer
- **IMPORTANT:** During installation, check the box "Add Python to PATH"
- Verify: Open Command Prompt (Windows) or Terminal (Mac) and type:
  ```
  python --version
  ```
  You should see: `Python 3.11.x` or higher

### Step 3: Open Project in VS Code

1. Open Visual Studio Code
2. Click **File → Open Folder** (or press `Ctrl+K Ctrl+O`)
3. Navigate to your `neurobotanica_project` folder
4. Click **Select Folder**

You should now see the project files in the left sidebar!

### Step 4: Install Python Extension

1. Click the **Extensions** icon on the left (looks like 4 squares)
2. Search for "Python"
3. Click **Install** on "Python" by Microsoft
4. Wait for installation to complete

### Step 5: Set Up Python Environment

1. In VS Code, open the **Terminal** (View → Terminal, or press `` Ctrl+` ``)
2. Type these commands one at a time:

**On Windows:**
```bash
python -m venv venv
venv\Scripts\activate
pip install --upgrade pip
pip install rdkit pandas numpy scikit-learn matplotlib seaborn pytest
```

**On Mac/Linux:**
```bash
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install rdkit pandas numpy scikit-learn matplotlib seaborn pytest
```

**Note:** Installation takes 5-10 minutes. RDKit is a large package.

You should see `(venv)` appear in your terminal prompt when activated.

### Step 6: Test Everything Works

In the terminal, run:
```bash
python src/analysis/cannabinoid_analyzer.py
```

You should see analysis output showing cannabinoid properties, receptor binding, and clinical studies!

### Step 7: Select Python Interpreter

1. Press `Ctrl+Shift+P` (Windows/Linux) or `Cmd+Shift+P` (Mac)
2. Type "Python: Select Interpreter"
3. Choose the one that shows `./venv/` or `venv`

Now VS Code will use your virtual environment for autocomplete and error checking!

---

## ✅ You're All Set!

**Your project is ready.** You can now:

- **View data:** Open `data/training/neurobotanica_complete_dataset_63compounds.json`
- **Edit code:** Open any `.py` file in `src/`
- **Run analysis:** `python src/analysis/cannabinoid_analyzer.py`
- **Run tests:** `pytest tests/ -v` (when you add tests)

---

## 🔧 Troubleshooting

**"python: command not found"**
- You need to install Python first (see Step 2B)
- Make sure you checked "Add Python to PATH" during installation

**"No module named 'rdkit'"**
- Your virtual environment isn't activated
- Look for `(venv)` at the start of your terminal line
- Activate it again:
  - Windows: `venv\Scripts\activate`
  - Mac/Linux: `source venv/bin/activate`

**"Cannot find file"**
- Make sure you're in the `neurobotanica_project` folder
- Check the terminal shows the right directory: `pwd` (Mac/Linux) or `cd` (Windows)

**RDKit won't install**
- Try using Anaconda instead: https://www.anaconda.com/download
- Then: `conda create -n neuro python=3.11`
- `conda activate neuro`
- `conda install -c conda-forge rdkit`
- `pip install pandas numpy scikit-learn matplotlib`

---

## 📚 Next Steps

1. **Read the full README:** Open `README.md` for detailed documentation
2. **Explore the data:** Look at JSON files in `data/training/`
3. **Modify code:** Edit scripts in `src/` to answer your research questions
4. **Add features:** Create new analysis scripts based on examples

---

## 💡 Daily Workflow

Each time you work on the project:

1. Open VS Code
2. Open your project folder (File → Open Folder)
3. Open Terminal (`` Ctrl+` ``)
4. Activate virtual environment:
   - Windows: `venv\Scripts\activate`
   - Mac/Linux: `source venv/bin/activate`
5. Start coding!

---

**Need help?** Check the full README.md for detailed documentation.
