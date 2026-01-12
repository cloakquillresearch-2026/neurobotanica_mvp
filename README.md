# NeuroBotanica - Cross-Kingdom Botanical Intelligence Platform

**Status:** Phase 1 Complete + Polysaccharide Expansion  
**Patent Filings:** Cannabis (Dec 22, 2025) + Polysaccharides (March 2026)  
**Version:** 1.1.0  
**Organization:** Cloak and Quill Research (501c3 Nonprofit)

---

## 🌿 Project Overview

NeuroBotanica is an AI-powered computational platform that predicts therapeutic compound synergies across **four botanical kingdoms** before physical testing:

### Botanical Kingdoms Supported
- **🌿 Cannabis Kingdom** - Dimeric cannabinoid prediction (e.g., cannabizetol)
- **🌊 Marine Polysaccharides** - Okra, kelp, seaweed therapeutic combinations
- **🍄 Fungal Polysaccharides** - Medicinal mushrooms, glucans, immunomodulatory compounds
- **🥕 Plant Polysaccharides** - Fenugreek, aloe, traditional botanical medicines

### Core Capabilities
- **Molecular Structure Prediction** - SAR modeling and therapeutic effect forecasting
- **Synergy Analysis** - Multi-compound optimization and entourage effects
- **Regulatory Compliance** - FDA-ready documentation generation
- **Traditional Knowledge Integration** - Ethical consent-gated traditional medicine validation
- **Computational Workflows** - Unified algorithms across all botanical domains

### Scientific Innovation
- **63+ cannabinoid compounds** with 40+ RDKit molecular descriptors each
- **41 clinical studies** from NORML database with condition-based inference
- **120+ dimeric cannabinoid predictions** with synergy scoring
- **Cross-kingdom polysaccharide analysis** with marine/fungal/plant applications
- **Patent-protected methodologies** (2 provisional filings)

### Business Model
- **SaaS Platform** - Subscription-based access for research institutions
- **API Services** - Enterprise integration for pharmaceutical companies
- **Nonprofit Revenue Distribution** - 70% to traditional knowledge communities
- **Market Opportunity** - $116B+ across four botanical kingdoms

---

## 📚 Documentation

### Core Documentation
- **[📖 Project Overview for Contributors](NeuroBotanica_Project_Overview.md)** - Comprehensive technical guide (640+ lines)
- **[🎓 Patent Teaching Guide](scripts/#%20📚%20NeuroBotanica%20Patent%20-%20Teaching%20Gui.md)** - Technical platform explanation
- **[🔬 Polysaccharide Expansion Strategy](scripts/#%20NeuroBotanica%20Polysaccharide%20Expansion.md)** - Cross-kingdom expansion details

### Patent Materials
- **[📄 Cannabis Provisional Patent](scripts/Provisional_Patent_NeuroBotanica_Polysaccharides.md)** - Filed December 22, 2025
- **[📄 Polysaccharide Provisional Patent](scripts/Provisional_Patent_NeuroBotanica_Polysaccharides.md)** - Target March 2026

### Deployment & Setup
- **[🚀 Quick Deployment Guide](DEPLOYMENT_QUICKSTART.md)** - Fast setup instructions
- **[☁️ Cloudflare Deployment](MANUAL_CLOUDFLARE_DEPLOYMENT.md)** - Edge computing deployment
- **[📋 Manual Deployment](MANUAL_DEPLOYMENT.md)** - Complete installation guide

### Additional Resources
- **[🏗️ MVP Overview](docs/NeuroBotanica_MVP_Overview.md)** - Product specifications
- **[📚 Deployment Guide](docs/NeuroBotanica_Deployment_Guide.md)** - Technical deployment details
- **[📊 Complete Update](docs/NeuroBotanica_Complete_Update.md)** - Project status report

---

## 🛠️ Getting Started

### Prerequisites

**Required Software:**
- **Python 3.11 or newer** (tested with 3.13) - Core runtime for backend and ML
- **Node.js 18+** - Frontend development and build tools
- **Git** - Version control system
- **VS Code** - Recommended IDE with extensions

**Recommended VS Code Extensions:**
- Python (Microsoft) - Language support and debugging
- Pylance (Microsoft) - Advanced Python IntelliSense
- TypeScript Importer - Auto-import for TypeScript
- Tailwind CSS IntelliSense - CSS class autocompletion
- GitLens - Enhanced Git capabilities
- Prettier - Code formatting
- ESLint - JavaScript/TypeScript linting

### Quick Setup

1. **Clone the Repository:**
```bash
git clone https://github.com/cloakquillresearch-2026/neurobotanica_mvp.git
cd neurobotanica_project
```

2. **Set Up Python Environment:**
```bash
# Create virtual environment
python -m venv venv

# Activate environment
# Windows
venv\Scripts\activate
# Mac/Linux
source venv/bin/activate

# Install dependencies
pip install --upgrade pip
pip install -r requirements.txt
```

3. **Set Up Node.js Environment:**
```bash
# Install frontend dependencies
cd frontend
npm install
# or
yarn install
```

4. **Configure Environment Variables:**
```bash
# Copy environment template
cp .env.example .env

# Edit .env with your configuration
# API keys, database URLs, etc.
```

5. **Run Development Servers:**
```bash
# Backend (from project root)
cd backend
uvicorn main:app --reload --host 0.0.0.0 --port 8000

# Frontend (from frontend directory)
npm run dev
# or
yarn dev
```

### Testing Setup

```bash
# Run all tests
pytest tests/ -v

# Run specific test files
pytest tests/test_installation.py
pytest tests/test_dimer_prediction.py
```

---

## 🏗️ Project Structure

```
neurobotanica_project/
├── backend/                 # FastAPI application
│   ├── api/                # API route handlers
│   ├── config/             # Configuration management
│   ├── middleware/         # Custom middleware (auth, CORS, etc.)
│   ├── models/             # Pydantic models and database schemas
│   ├── routers/            # API route definitions
│   ├── services/           # Business logic and external integrations
│   │   ├── biopath/       # Biological pathway analysis
│   │   ├── chempath/      # Chemical pathway modeling
│   │   └── toxpath/       # Toxicity assessment
│   ├── tests/              # Backend unit and integration tests
│   └── utils/              # Utility functions and helpers
├── frontend/               # Next.js application
│   ├── app/               # App Router pages and layouts
│   ├── components/        # Reusable React components
│   ├── lib/               # Utility functions and configurations
│   ├── public/            # Static assets
│   └── styles/            # Global styles and Tailwind config
├── data/                  # Data processing and storage
│   ├── clinical_evidence/ # Clinical trial data
│   ├── expanded_studies/  # Processed research datasets
│   ├── norml_extraction/  # Cannabis research database
│   └── processed/         # Cleaned and formatted datasets
├── models/                # Machine learning models and checkpoints
├── tests/                 # Integration and end-to-end tests
├── scripts/               # Utility scripts and automation
│   ├── *.py              # Python automation scripts
│   ├── *.bat             # Windows batch scripts
│   └── *.sh              # Shell scripts
├── docs/                  # Documentation and guides
└── runs/                  # Experiment logs and results
```

---

## 🤝 Contributing

We welcome contributions from developers, researchers, and domain experts! Here's how to get involved:

### Ways to Contribute

1. **Code Development** - Backend APIs, frontend components, ML pipelines
2. **Documentation** - Improve guides, add tutorials, create examples
3. **Testing** - Write unit tests, integration tests, performance benchmarks
4. **Research** - Validate algorithms, analyze results, suggest improvements
5. **Design** - UI/UX improvements, data visualization, user experience

### Getting Started with Contributions

1. **Read the Contributor Guide** - See [Project Overview](NeuroBotanica_Project_Overview.md) for detailed technical tasks
2. **Check Open Issues** - Review GitHub issues for current priorities
3. **Start Small** - Begin with documentation or simple bug fixes
4. **Follow Guidelines** - Use conventional commits, write tests, document changes

### Development Workflow

```bash
# 1. Create feature branch
git checkout -b feature/your-feature-name

# 2. Make changes and test
# ... development work ...

# 3. Run tests
pytest backend/tests/
npm test

# 4. Commit changes
git add .
git commit -m "feat: add your feature description"

# 5. Push and create PR
git push origin feature/your-feature-name
```

### Code Quality Standards

- **Python:** PEP 8 compliance, type hints, comprehensive docstrings
- **JavaScript/TypeScript:** ESLint rules, Prettier formatting
- **Testing:** 80%+ code coverage, integration tests for APIs
- **Documentation:** All public APIs documented, inline comments for complex logic

---

## 🗺️ Project Roadmap

### Phase 1 (Completed)
- ✅ Core dataset and ML framework
- ✅ Cannabis dimeric prediction algorithms
- ✅ Basic web platform and API
- ✅ Provisional patent filing (Cannabis)

### Phase 2 (Current - Q1 2026)
- 🔄 Polysaccharide expansion integration
- 🔄 Advanced visualization features
- 🔄 Beta API launch for research partners
- 🔄 Provisional patent filing (Polysaccharides)

### Phase 3 (2026)
- 📅 Full SaaS platform launch
- 📅 Mobile application development
- 📅 Integration with major cheminformatics databases
- 📅 Multi-language support (R, Julia bindings)

### Phase 4 (2027+)
- 📅 Global research network integration
- 📅 Real-time collaborative analysis tools
- 📅 Advanced AI model development
- 📅 Expansion to additional therapeutic domains

---

## 🧪 Technical Architecture

### Backend Stack
- **Framework:** FastAPI (async Python web framework)
- **Database:** Cloudflare D1 + KV storage
- **Authentication:** Firebase Auth with RBAC
- **Deployment:** Cloudflare Workers (<200ms global response times)

### Frontend Stack
- **Framework:** Next.js 14 with App Router
- **Styling:** Tailwind CSS + Radix UI components
- **State Management:** React hooks + Context API
- **Visualization:** D3.js + Recharts for scientific data

### Machine Learning Stack
- **Libraries:** scikit-learn, RDKit, XGBoost
- **Models:** SAR models, ensemble methods, predictive algorithms
- **Data Processing:** Pandas, NumPy, SciPy
- **Validation:** Cross-validation, hyperparameter tuning

### DevOps & Infrastructure
- **Version Control:** Git + GitHub
- **CI/CD:** GitHub Actions
- **Testing:** pytest + coverage reporting
- **Monitoring:** Cloudflare Analytics + custom logging

---

## 📊 Performance Benchmarks

- **API Response Time:** <100ms for simple queries
- **Model Prediction Time:** <5 seconds for molecular analysis
- **Concurrent Users:** 1000+ supported
- **Data Processing:** 10,000+ molecules per hour
- **Test Coverage:** Target 90%+

---

## 📞 Support & Community

### Communication Channels
- **GitHub Issues:** Bug reports and feature requests
- **Pull Request Discussions:** Code review and collaboration
- **Discord/Slack:** Real-time community discussions
- **Monthly Calls:** Project updates and community meetings

### Getting Help
- **Documentation:** Comprehensive guides in `/docs` and `/scripts`
- **Office Hours:** Regular availability for technical questions
- **Mentorship:** Pair programming and code review support
- **Tech Talks:** Weekly sessions on relevant technologies

### Recognition
- Contributor acknowledgments in release notes
- Co-authorship opportunities on publications
- Speaking opportunities at conferences
- Priority access to new features

---

## 📜 License & Legal

**Organization:** Cloak and Quill Research (501c3 Nonprofit)  
**License:** Proprietary (research use only)  
**Patents:** Provisional filings for core methodologies  
**Compliance:** HIPAA, FDA regulatory frameworks  

**Important:** This codebase contains trade secrets and patent-pending technologies. Access is restricted to authorized contributors only.

---

## 🙏 Acknowledgments

- **Scientific Advisors:** Domain experts in cheminformatics and pharmacology
- **Open Source Community:** Contributors to RDKit, scikit-learn, and related libraries
- **Research Partners:** Academic institutions and pharmaceutical companies
- **Traditional Knowledge Communities:** Indigenous and traditional medicine practitioners

---

## 🚀 Quick Start Commands

```bash
# Clone and setup
git clone https://github.com/cloakquillresearch-2026/neurobotanica_mvp.git
cd neurobotanica_project

# Setup environment
python -m venv venv
venv\Scripts\activate  # Windows
pip install -r requirements.txt

# Start development
cd backend && uvicorn main:app --reload --host 0.0.0.0 --port 8000
# In another terminal: cd frontend && npm run dev

# Run tests
pytest backend/tests/
```

---

**Ready to contribute?** Start with our [Project Overview](NeuroBotanica_Project_Overview.md) and join us in building the future of computational botanical therapeutics! 🌿🚀
