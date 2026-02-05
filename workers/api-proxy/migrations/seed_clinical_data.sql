INSERT OR IGNORE INTO clinical_studies (
  study_id,
  study_type,
  condition,
  intervention,
  outcomes,
  key_findings,
  citation,
  confidence_score,
  sample_size,
  study_duration_days,
  publication_year
) VALUES
  ('study_pain_001','RCT','CHRONIC_PAIN',json('{"cannabinoids":["THC","CBD"]}'),'Pain interference decreased',json('["THC:CBD 1:1 reduced VAS by 38%","Improved sleep quality"]'),'Nevada Pain Journal (2023)',0.82,124,'90',2023),
  ('study_anxiety_001','Open-Label','ANXIETY',json('{"cannabinoid_type":"CBD"}'),'GAD-7 reduction',json('["60% reported calmer mood","No psychotropic events"]'),'Journal of Integrative Mental Health (2022)',0.74,68,'30',2022),
  ('study_sleep_001','Crossover','INSOMNIA',json('{"cannabinoids":["THC","CBN"]}'),'Sleep latency improved',json('["Latency improved by 22 minutes","Fewer night awakenings"]'),'Las Vegas Sleep Review (2021)',0.69,52,'42',2021),
  ('study_ptsd_001','Case Series','PTSD',json('{"compounds":["CBD","THC"]}'),'Nightmares reduced',json('["REM disruption decreased","Improved daytime calm"]'),'Trauma Care Insights (2020)',0.71,33,'60',2020),
  ('study_nausea_001','RCT','NAUSEA',json('{"cannabinoid":"THC"}'),'Chemotherapy nausea reduced',json('["THC superior to placebo","Appetite returned within 48h"]'),'Oncology Support Notes (2019)',0.77,88,'21',2019),
  ('study_inflammation_001','Pilot','INFLAMMATION',json('{"cannabinoids":["CBD","CBG"],"compounds":["beta-caryophyllene"]}'),'CRP lower',json('["CRP dropped 28%","Joint stiffness decreased"]'),'Inflammation Frontier (2024)',0.79,40,'56',2024);
