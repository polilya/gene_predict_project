from dataclasses import dataclass

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
    final_data: str

@dataclass
class Params:
    num_of_samples: int

@dataclass
class MLConfig:
    paths: Paths
    files: Files
    params: Params