"""
ML Models - NeuroBotanica Week 6
Therapeutic prediction models with uncertainty quantification.

Models:
1. TherapeuticPredictionModel - Condition → Optimal cannabinoid efficacy
2. PatientResponseModel - Patient → Treatment response probability  
3. DimerPotentialModel - Dimer → Therapeutic potential

All models include:
- Confidence-weighted training
- Uncertainty quantification
- Cross-validation metrics
- Model persistence
"""
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
import logging
import json

# ML imports
from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier, RandomForestRegressor
from sklearn.model_selection import cross_val_score, train_test_split, KFold
from sklearn.metrics import r2_score, mean_squared_error, accuracy_score, roc_auc_score, classification_report
from sklearn.preprocessing import StandardScaler, LabelEncoder
import joblib

logger = logging.getLogger(__name__)


@dataclass
class ModelMetrics:
    """Training and evaluation metrics for a model."""
    model_name: str
    model_version: str
    train_score: float
    test_score: float
    cv_scores: List[float]
    cv_mean: float
    cv_std: float
    num_features: int
    num_samples: int
    training_time_seconds: float
    feature_importances: Dict[str, float]
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())
    
    def to_dict(self) -> Dict:
        return {
            "model_name": self.model_name,
            "model_version": self.model_version,
            "train_score": round(self.train_score, 4),
            "test_score": round(self.test_score, 4),
            "cv_mean": round(self.cv_mean, 4),
            "cv_std": round(self.cv_std, 4),
            "num_features": self.num_features,
            "num_samples": self.num_samples,
            "training_time_seconds": round(self.training_time_seconds, 2),
            "top_features": dict(sorted(self.feature_importances.items(), key=lambda x: x[1], reverse=True)[:10]),
            "created_at": self.created_at
        }


@dataclass
class PredictionResult:
    """Result from model prediction with uncertainty."""
    prediction: Union[float, np.ndarray]
    uncertainty: Union[float, np.ndarray]
    confidence_interval_lower: Union[float, np.ndarray]
    confidence_interval_upper: Union[float, np.ndarray]
    confidence_level: float = 0.95
    
    def to_dict(self) -> Dict:
        pred = self.prediction
        if isinstance(pred, np.ndarray):
            pred = pred.tolist() if len(pred) > 1 else float(pred[0])
        
        unc = self.uncertainty
        if isinstance(unc, np.ndarray):
            unc = unc.tolist() if len(unc) > 1 else float(unc[0])
        
        ci_lower = self.confidence_interval_lower
        if isinstance(ci_lower, np.ndarray):
            ci_lower = ci_lower.tolist() if len(ci_lower) > 1 else float(ci_lower[0])
        
        ci_upper = self.confidence_interval_upper
        if isinstance(ci_upper, np.ndarray):
            ci_upper = ci_upper.tolist() if len(ci_upper) > 1 else float(ci_upper[0])
        
        return {
            "prediction": pred,
            "uncertainty": unc,
            "confidence_interval": [ci_lower, ci_upper],
            "confidence_level": self.confidence_level
        }


class BaseModel:
    """Base class for NeuroBotanica ML models."""
    
    MODEL_VERSION = "1.0.0"
    
    def __init__(self, model_name: str):
        self.model_name = model_name
        self.model = None
        self.scaler = StandardScaler()
        self.feature_names: List[str] = []
        self.is_trained = False
        self.metrics: Optional[ModelMetrics] = None
    
    def _get_feature_importances(self) -> Dict[str, float]:
        """Get feature importance scores."""
        if not hasattr(self.model, 'feature_importances_'):
            return {}
        
        importances = self.model.feature_importances_
        return {
            name: float(imp) 
            for name, imp in zip(self.feature_names, importances)
        }
    
    def save(self, path: str) -> None:
        """Save model to disk."""
        save_data = {
            "model": self.model,
            "scaler": self.scaler,
            "feature_names": self.feature_names,
            "metrics": self.metrics.to_dict() if self.metrics else None,
            "model_name": self.model_name,
            "model_version": self.MODEL_VERSION,
            "saved_at": datetime.now().isoformat()
        }
        
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(save_data, path)
        logger.info(f"Model saved to {path}")
    
    def load(self, path: str) -> None:
        """Load model from disk."""
        save_data = joblib.load(path)
        
        self.model = save_data["model"]
        self.scaler = save_data["scaler"]
        self.feature_names = save_data["feature_names"]
        self.is_trained = True
        
        logger.info(f"Model loaded from {path}")


class TherapeuticPredictionModel(BaseModel):
    """Gradient boosting model for therapeutic efficacy prediction.
    
    Predicts efficacy score (0-1) for cannabinoid-condition combinations.
    Uses 2D/3D molecular descriptors and receptor affinities.
    """
    
    def __init__(
        self,
        n_estimators: int = 100,
        learning_rate: float = 0.1,
        max_depth: int = 5,
        random_state: int = 42
    ):
        super().__init__("TherapeuticPrediction")
        
        self.model = GradientBoostingRegressor(
            n_estimators=n_estimators,
            learning_rate=learning_rate,
            max_depth=max_depth,
            random_state=random_state,
            subsample=0.7,              # Increased regularization (was 0.8)
            min_samples_split=10,       # Increased from 5 to prevent overfitting
            min_samples_leaf=5,         # Increased from 2 to require more samples per leaf
            max_features='sqrt'         # Add feature subsampling for regularization
        )
        self.random_state = random_state
    
    def train(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        sample_weights: Optional[pd.Series] = None,
        test_size: float = 0.2
    ) -> ModelMetrics:
        """Train model with optional confidence weights.
        
        Args:
            X: Feature DataFrame
            y: Target values (efficacy scores 0-1)
            sample_weights: Optional confidence weights
            test_size: Fraction of data for testing
            
        Returns:
            ModelMetrics with training results
        """
        import time
        start_time = time.time()
        
        self.feature_names = list(X.columns)
        
        # Scale features
        X_scaled = self.scaler.fit_transform(X)
        
        # Split data
        if sample_weights is not None:
            X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(
                X_scaled, y, sample_weights, 
                test_size=test_size, 
                random_state=self.random_state
            )
        else:
            X_train, X_test, y_train, y_test = train_test_split(
                X_scaled, y, 
                test_size=test_size, 
                random_state=self.random_state
            )
            w_train, w_test = None, None
        
        # Train model
        self.model.fit(X_train, y_train, sample_weight=w_train)
        
        # Evaluate
        train_score = self.model.score(X_train, y_train)
        test_score = self.model.score(X_test, y_test)
        
        # Cross-validation
        cv_scores = cross_val_score(
            self.model, X_scaled, y, 
            cv=5, scoring='r2'
        )
        
        training_time = time.time() - start_time
        
        self.metrics = ModelMetrics(
            model_name=self.model_name,
            model_version=self.MODEL_VERSION,
            train_score=train_score,
            test_score=test_score,
            cv_scores=cv_scores.tolist(),
            cv_mean=cv_scores.mean(),
            cv_std=cv_scores.std(),
            num_features=len(self.feature_names),
            num_samples=len(X),
            training_time_seconds=training_time,
            feature_importances=self._get_feature_importances()
        )
        
        self.is_trained = True
        logger.info(f"Model trained: R²={test_score:.4f}, CV={cv_scores.mean():.4f}±{cv_scores.std():.4f}")
        
        return self.metrics
    
    def predict(self, X: pd.DataFrame) -> np.ndarray:
        """Simple prediction without uncertainty."""
        if not self.is_trained:
            raise ValueError("Model not trained")
        
        X_scaled = self.scaler.transform(X)
        return self.model.predict(X_scaled)
    
    def predict_with_uncertainty(self, X: pd.DataFrame) -> PredictionResult:
        """Predict with uncertainty quantification using ensemble variance.
        
        Uses individual estimator predictions to estimate uncertainty.
        """
        if not self.is_trained:
            raise ValueError("Model not trained")
        
        X_scaled = self.scaler.transform(X)
        
        # Get predictions from each estimator (tree)
        estimator_predictions = np.array([
            estimator[0].predict(X_scaled)
            for estimator in self.model.estimators_
        ])
        
        # Calculate mean and std
        mean_pred = estimator_predictions.mean(axis=0)
        std_pred = estimator_predictions.std(axis=0)
        
        # Clip predictions to valid range
        mean_pred = np.clip(mean_pred, 0.0, 1.0)
        
        # 95% confidence interval
        z_score = 1.96
        ci_lower = np.clip(mean_pred - z_score * std_pred, 0.0, 1.0)
        ci_upper = np.clip(mean_pred + z_score * std_pred, 0.0, 1.0)
        
        return PredictionResult(
            prediction=mean_pred,
            uncertainty=std_pred,
            confidence_interval_lower=ci_lower,
            confidence_interval_upper=ci_upper
        )
    
    def get_condition_recommendations(
        self,
        X: pd.DataFrame,
        conditions: List[str],
        top_n: int = 5
    ) -> List[Dict]:
        """Get top condition recommendations for given features."""
        result = self.predict_with_uncertainty(X)
        
        recommendations = []
        for i, condition in enumerate(conditions):
            if i < len(result.prediction):
                recommendations.append({
                    "condition": condition,
                    "efficacy": float(result.prediction[i]),
                    "uncertainty": float(result.uncertainty[i]),
                    "confidence_interval": [
                        float(result.confidence_interval_lower[i]),
                        float(result.confidence_interval_upper[i])
                    ]
                })
        
        # Sort by efficacy
        recommendations.sort(key=lambda x: x["efficacy"], reverse=True)
        return recommendations[:top_n]


class PatientResponseModel(BaseModel):
    """Classification model for patient treatment response prediction.
    
    Predicts probability of positive response (binary classification).
    """
    
    def __init__(
        self,
        n_estimators: int = 100,
        learning_rate: float = 0.1,
        max_depth: int = 4,
        random_state: int = 42
    ):
        super().__init__("PatientResponse")
        
        self.model = GradientBoostingClassifier(
            n_estimators=n_estimators,
            learning_rate=learning_rate,
            max_depth=max_depth,
            random_state=random_state,
            subsample=0.8
        )
        self.random_state = random_state
    
    def train(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        sample_weights: Optional[pd.Series] = None,
        test_size: float = 0.2
    ) -> ModelMetrics:
        """Train binary classification model."""
        import time
        start_time = time.time()
        
        self.feature_names = list(X.columns)
        
        # Scale features
        X_scaled = self.scaler.fit_transform(X)
        
        # Ensure binary target
        y_binary = (y > 0.5).astype(int)
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X_scaled, y_binary,
            test_size=test_size,
            random_state=self.random_state,
            stratify=y_binary
        )
        
        # Train
        self.model.fit(X_train, y_train)
        
        # Evaluate
        train_pred = self.model.predict(X_train)
        test_pred = self.model.predict(X_test)
        
        train_score = accuracy_score(y_train, train_pred)
        test_score = accuracy_score(y_test, test_pred)
        
        # Cross-validation
        cv_scores = cross_val_score(
            self.model, X_scaled, y_binary,
            cv=5, scoring='accuracy'
        )
        
        training_time = time.time() - start_time
        
        self.metrics = ModelMetrics(
            model_name=self.model_name,
            model_version=self.MODEL_VERSION,
            train_score=train_score,
            test_score=test_score,
            cv_scores=cv_scores.tolist(),
            cv_mean=cv_scores.mean(),
            cv_std=cv_scores.std(),
            num_features=len(self.feature_names),
            num_samples=len(X),
            training_time_seconds=training_time,
            feature_importances=self._get_feature_importances()
        )
        
        self.is_trained = True
        logger.info(f"Model trained: Accuracy={test_score:.4f}, CV={cv_scores.mean():.4f}")
        
        return self.metrics
    
    def predict_proba(self, X: pd.DataFrame) -> np.ndarray:
        """Predict response probabilities."""
        if not self.is_trained:
            raise ValueError("Model not trained")
        
        X_scaled = self.scaler.transform(X)
        return self.model.predict_proba(X_scaled)[:, 1]
    
    def predict_with_uncertainty(self, X: pd.DataFrame) -> PredictionResult:
        """Predict with uncertainty using staged predictions."""
        if not self.is_trained:
            raise ValueError("Model not trained")
        
        X_scaled = self.scaler.transform(X)
        
        # Get staged predictions for uncertainty
        staged_probs = np.array([
            pred[:, 1] for pred in self.model.staged_predict_proba(X_scaled)
        ])
        
        # Use last N stages for uncertainty estimation
        n_stages = min(20, len(staged_probs))
        recent_staged = staged_probs[-n_stages:]
        
        mean_prob = recent_staged.mean(axis=0)
        std_prob = recent_staged.std(axis=0)
        
        return PredictionResult(
            prediction=mean_prob,
            uncertainty=std_prob,
            confidence_interval_lower=np.clip(mean_prob - 1.96 * std_prob, 0, 1),
            confidence_interval_upper=np.clip(mean_prob + 1.96 * std_prob, 0, 1)
        )


class DimerPotentialModel(BaseModel):
    """Model for predicting dimer therapeutic potential.
    
    Uses combined features from both parent compounds plus dimer-specific features.
    """
    
    def __init__(
        self,
        n_estimators: int = 150,
        learning_rate: float = 0.08,
        max_depth: int = 6,
        random_state: int = 42
    ):
        super().__init__("DimerPotential")
        
        self.model = GradientBoostingRegressor(
            n_estimators=n_estimators,
            learning_rate=learning_rate,
            max_depth=max_depth,
            random_state=random_state,
            subsample=0.8,
            min_samples_split=3
        )
        self.random_state = random_state
    
    def train(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        sample_weights: pd.Series = None,
        test_size: float = 0.2
    ) -> ModelMetrics:
        """Train dimer potential prediction model."""
        import time
        start_time = time.time()
        self.feature_names = list(X.columns)
        # Scale features
        X_scaled = self.scaler.fit_transform(X)
        # Split
        X_train, X_test, y_train, y_test = train_test_split(
            X_scaled, y,
            test_size=test_size,
            random_state=self.random_state
        )
        # Train
        if sample_weights is not None:
            w_train = sample_weights.iloc[X_train.index] if hasattr(X_train, 'index') else sample_weights[:len(X_train)]
            self.model.fit(X_train, y_train, sample_weight=w_train)
        else:
            self.model.fit(X_train, y_train)
        # Evaluate
        train_score = self.model.score(X_train, y_train)
        test_score = self.model.score(X_test, y_test)
        # Cross-validation
        from sklearn.model_selection import KFold
        if sample_weights is not None:
            # Manual cross-validation to support sample weights
            kf = KFold(n_splits=5, shuffle=True, random_state=self.random_state)
            cv_scores = []
            for train_idx, test_idx in kf.split(X_scaled):
                X_tr, X_te = X_scaled[train_idx], X_scaled[test_idx]
                y_tr, y_te = y.iloc[train_idx], y.iloc[test_idx]
                w_tr = sample_weights.iloc[train_idx]
                model_cv = GradientBoostingRegressor(
                    n_estimators=self.model.n_estimators,
                    learning_rate=self.model.learning_rate,
                    max_depth=self.model.max_depth,
                    random_state=self.random_state,
                    subsample=self.model.subsample,
                    min_samples_split=self.model.min_samples_split
                )
                model_cv.fit(X_tr, y_tr, sample_weight=w_tr)
                score = model_cv.score(X_te, y_te)
                cv_scores.append(score)
            cv_scores = np.array(cv_scores)
        else:
            cv_scores = cross_val_score(
                self.model, X_scaled, y,
                cv=5, scoring='r2'
            )
        training_time = time.time() - start_time
        self.metrics = ModelMetrics(
            model_name=self.model_name,
            model_version=self.MODEL_VERSION,
            train_score=train_score,
            test_score=test_score,
            cv_scores=cv_scores.tolist(),
            cv_mean=cv_scores.mean(),
            cv_std=cv_scores.std(),
            num_features=len(self.feature_names),
            num_samples=len(X),
            training_time_seconds=training_time,
            feature_importances=self._get_feature_importances()
        )
        self.is_trained = True
        logger.info(f"Dimer model trained: R²={test_score:.4f}")
        return self.metrics
    
    def predict_with_uncertainty(self, X: pd.DataFrame) -> PredictionResult:
        """Predict dimer potential with uncertainty."""
        if not self.is_trained:
            raise ValueError("Model not trained")
        
        X_scaled = self.scaler.transform(X)
        
        # Ensemble predictions
        estimator_predictions = np.array([
            estimator[0].predict(X_scaled)
            for estimator in self.model.estimators_
        ])
        
        mean_pred = estimator_predictions.mean(axis=0)
        std_pred = estimator_predictions.std(axis=0)
        
        mean_pred = np.clip(mean_pred, 0.0, 1.0)
        
        return PredictionResult(
            prediction=mean_pred,
            uncertainty=std_pred,
            confidence_interval_lower=np.clip(mean_pred - 1.96 * std_pred, 0, 1),
            confidence_interval_upper=np.clip(mean_pred + 1.96 * std_pred, 0, 1)
        )
    
    def rank_dimers(
        self,
        dimers: List[Dict],
        X: pd.DataFrame
    ) -> List[Dict]:
        """Rank dimers by predicted therapeutic potential."""
        result = self.predict_with_uncertainty(X)
        
        ranked = []
        for i, dimer in enumerate(dimers):
            if i < len(result.prediction):
                ranked.append({
                    "dimer_name": dimer.get("dimer_name", f"Dimer_{i}"),
                    "therapeutic_potential": float(result.prediction[i]),
                    "uncertainty": float(result.uncertainty[i]),
                    "confidence_interval": [
                        float(result.confidence_interval_lower[i]),
                        float(result.confidence_interval_upper[i])
                    ],
                    "parent_1": dimer.get("parent_1_name", ""),
                    "parent_2": dimer.get("parent_2_name", "")
                })
        
        ranked.sort(key=lambda x: x["therapeutic_potential"], reverse=True)
        return ranked


class ModelRegistry:
    """Registry for managing trained models."""
    
    def __init__(self, models_dir: str = "models"):
        self.models_dir = Path(models_dir)
        self.models_dir.mkdir(parents=True, exist_ok=True)
        self.registry_file = self.models_dir / "registry.json"
        self.registry: Dict[str, Dict] = self._load_registry()
    
    def _load_registry(self) -> Dict:
        """Load registry from disk."""
        if self.registry_file.exists():
            with open(self.registry_file, 'r') as f:
                return json.load(f)
        return {}
    
    def _save_registry(self) -> None:
        """Save registry to disk."""
        with open(self.registry_file, 'w') as f:
            json.dump(self.registry, f, indent=2)
    
    def register_model(
        self,
        model: BaseModel,
        tags: Optional[List[str]] = None
    ) -> str:
        """Register and save a trained model."""
        if not model.is_trained:
            raise ValueError("Cannot register untrained model")
        
        # Generate unique ID
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        model_id = f"{model.model_name}_v{model.MODEL_VERSION}_{timestamp}"
        
        # Save model
        model_path = self.models_dir / f"{model_id}.joblib"
        model.save(str(model_path))
        
        # Register metadata
        self.registry[model_id] = {
            "model_name": model.model_name,
            "model_version": model.MODEL_VERSION,
            "path": str(model_path),
            "metrics": model.metrics.to_dict() if model.metrics else None,
            "tags": tags or [],
            "registered_at": datetime.now().isoformat()
        }
        
        self._save_registry()
        logger.info(f"Model registered: {model_id}")
        
        return model_id
    
    def get_model(self, model_id: str) -> BaseModel:
        """Load a registered model."""
        if model_id not in self.registry:
            raise ValueError(f"Model not found: {model_id}")
        
        model_info = self.registry[model_id]
        model_name = model_info["model_name"]
        
        # Create appropriate model instance
        if model_name == "TherapeuticPrediction":
            model = TherapeuticPredictionModel()
        elif model_name == "PatientResponse":
            model = PatientResponseModel()
        elif model_name == "DimerPotential":
            model = DimerPotentialModel()
        else:
            raise ValueError(f"Unknown model type: {model_name}")
        
        model.load(model_info["path"])
        return model
    
    def list_models(
        self,
        model_name: Optional[str] = None,
        tags: Optional[List[str]] = None
    ) -> List[Dict]:
        """List registered models with optional filtering."""
        results = []
        
        for model_id, info in self.registry.items():
            # Filter by name
            if model_name and info["model_name"] != model_name:
                continue
            
            # Filter by tags
            if tags:
                if not any(t in info.get("tags", []) for t in tags):
                    continue
            
            results.append({
                "model_id": model_id,
                **info
            })
        
        return sorted(results, key=lambda x: x["registered_at"], reverse=True)
    
    def get_best_model(self, model_name: str, metric: str = "test_score") -> Optional[str]:
        """Get best model ID by metric."""
        models = self.list_models(model_name=model_name)
        
        if not models:
            return None
        
        best = max(
            models,
            key=lambda x: x.get("metrics", {}).get(metric, 0)
        )
        
        return best["model_id"]
