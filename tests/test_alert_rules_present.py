import os
import yaml


def test_prometheus_rules_exist_and_parsable():
    path = os.path.join('monitoring', 'prometheus', 'rules', 'inflammatory_rules.yml')
    assert os.path.exists(path), 'Prometheus rules file is missing'
    with open(path) as fh:
        data = yaml.safe_load(fh)
    assert 'groups' in data
    found = False
    for g in data.get('groups', []):
        for r in g.get('rules', []):
            if r.get('alert') == 'NeuroBotanicaInflammatoryHighLatency':
                found = True
    assert found, 'Expected alert NeuroBotanicaInflammatoryHighLatency not found'