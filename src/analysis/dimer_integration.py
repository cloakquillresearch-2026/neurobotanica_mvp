from typing import List, Union
from .dimer_schema import DimerEvidence

def merge_dimeric_evidence(training_data: List[Union[dict, DimerEvidence]], dimeric_data: List[Union[dict, DimerEvidence]]) -> List[Union[dict, DimerEvidence]]:
    merged = training_data.copy()
    for dimer in dimeric_data:
        if dimer not in merged:
            merged.append(dimer)
    return merged
