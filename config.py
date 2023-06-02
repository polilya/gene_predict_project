from dataclasses import dataclass

@dataclass
class ModelParams:
    class_weight: dict
    kernel: str
    gamma: str
    decision_function_shape: str
    C: float
    tol: float
    probability: bool

@dataclass
class AddParams:
    windows: list
    step: int
    plot: bool
    
@dataclass
class PPParams:
    num_of_samples: int
    start_min_length: int
    start_max_length: int
    end_min_length: int
    end_max_length: int

@dataclass
class Paths:
    log: str
    data: str

@dataclass
class Files:
    sequence: str
    cds: str
    cds_location: str
    exons: str
    features: str
    first_cds_location: str
    last_cds_location: str
    start_features: str
    end_features: str

@dataclass
class MLConfig:
    paths: Paths
    files: Files
    params: ModelParams
    add_params: AddParams
    preprocess_params: PPParams

