"""Targeted Google Scholar queries via SerpApi, 2022-2026. Logs every query and result to scholar_results.json."""
import json, os, sys, time
from pathlib import Path
import requests

HERE = Path(__file__).parent
QUERIES = [
    'online environment stiffness estimation robot force control',
    'online contact stiffness identification time-varying environment robot',
    'environment stiffness estimation Kalman filter robot interaction unknown surface',
    'contact parameter estimation moving environment robot force control',
    'soft tissue stiffness online estimation robot force position',
    'robotic ultrasound tissue stiffness estimation force control',
    'Hunt-Crossley online identification robot soft tissue',
    'viscoelastic parameter online estimation robot palpation extended Kalman filter',
    'environment stiffness estimation physiological motion compensation robot',
    'stiffness estimation excitation signal robot contact force recursive least squares',
    'interaction force and environment stiffness joint estimation robot',
    'adaptive force control unknown time-varying stiffness estimation robot manipulator',
    'robust stiffness estimation outliers contact robot',
    'moving horizon estimation environment stiffness robot',
    'Gaussian process contact stiffness estimation robot online',
    'learning-based environment stiffness estimation robot contact neural network',
]


def run(q, start=0):
    r = requests.get('https://serpapi.com/search.json', params=dict(engine='google_scholar', q=q, as_ylo=2022, num=20, start=start,
                                                                     api_key=os.environ['SERPAPI_KEY']), timeout=60)
    r.raise_for_status()
    return r.json().get('organic_results', [])


def main():
    out = []
    for q in QUERIES:
        try:
            res = run(q)
        except Exception as e:
            out.append(dict(query=q, error=str(e))); continue
        out.append(dict(query=q, results=[dict(title=x.get('title'), link=x.get('link'), info=(x.get('publication_info') or {}).get('summary'),
                                              cited=((x.get('inline_links') or {}).get('cited_by') or {}).get('total'), snippet=x.get('snippet'))
                                         for x in res]))
        time.sleep(1)
    (HERE/'scholar_results.json').write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding='utf-8')
    for o in out:
        print('##', o['query'], o.get('error', ''))
        for x in o.get('results', []):
            print('  -', (x['title'] or '')[:110], '|', (x['info'] or '')[:90], '| cited', x['cited'])


if __name__ == '__main__':
    main()
