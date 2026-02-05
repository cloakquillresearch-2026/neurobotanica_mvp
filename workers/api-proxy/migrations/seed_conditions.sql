INSERT OR IGNORE INTO conditions (
  condition_id,
  condition_name,
  category,
  recommended_cannabinoids,
  evidence_count,
  description,
  updated_at
) VALUES
  ('chronic_pain','Chronic Pain','pain','THC,CBD',25,'Evidence-backed analgesic guidance',datetime('now')),
  ('anxiety','Anxiety','anxiety','CBD,CBG',18,'Calming protocols with beta-caryophyllene support',datetime('now')),
  ('insomnia','Insomnia','sleep','THC,CBN',14,'Sleep-onset and maintenance optimization',datetime('now')),
  ('ptsd','PTSD','trauma','CBD,THC',12,'Trauma-informed cannabinoid protocols',datetime('now')),
  ('nausea','Nausea','general','THC,THCV',9,'Antiemetic cannabinoid support profiles',datetime('now')),
  ('inflammation','Inflammation','inflammation','CBD,CBG',22,'Anti-inflammatory composite guidance',datetime('now'));
