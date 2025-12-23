import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LogisticRegression, Ridge
from src.analysis.nextgen_dimer_ingest import load_nextgen_dimer_evidence
from src.ml_models.nextgen_dimer_input import encode_dimer_features
from src.ml_models.nextgen_dimer_fusion import prepare_batch

if __name__ == "__main__":
    # Load data
    dimers = load_nextgen_dimer_evidence("data/processed/nextgen_dimers.json")
    entries = [d.__dict__ for d in dimers]
    X = prepare_batch(entries)
    y_synergy = np.array([1 if e['interaction_effect'] == 'synergistic' else 0 for e in entries])
    y_effect = np.array([float(e.get('synergy_score', 1.0)) for e in entries])

    # Hyperparameter grid for LogisticRegression (synergy)
    synergy_grid = {
        'C': [0.01, 0.1, 1, 10],
        'solver': ['liblinear', 'lbfgs']
    }
    synergy_clf = GridSearchCV(LogisticRegression(), synergy_grid, cv=2)
    synergy_clf.fit(X, y_synergy)

    # Hyperparameter grid for Ridge (effect size)
    effect_grid = {
        'alpha': [0.01, 0.1, 1, 10, 100]
    }
    effect_reg = GridSearchCV(Ridge(), effect_grid, cv=2)
    effect_reg.fit(X, y_effect)

    print("Best synergy classifier params:", synergy_clf.best_params_)
    print("Best effect regressor params:", effect_reg.best_params_)
    print("Best synergy score (mean CV):", synergy_clf.best_score_)
    print("Best effect score (mean CV):", effect_reg.best_score_)
