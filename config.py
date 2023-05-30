from dataclasses import dataclass

@dataclass
class ModelParams:
    class_weight: int
    windows: list
    left_boundary: int
    right_boundary: int
    step: int
    
@dataclass
class PPParams:
    num_of_samples: int

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

@dataclass
class MLConfig:
    paths: Paths
    files: Files
    model_params: ModelParams
    preprocess_params: PPParams