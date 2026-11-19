"""Audit every alpha/beta force factor against its declared analytic formula."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
import jsbsim

COEFF_SINGLE=25.293749861383475
COEFF_STACK=52.00463519814777
KINDS=('Drag_crossflow','Drag_beta','Lift_alpha','Side_beta','Lift_crossflow','Side_crossflow')


def patch(root):
    for aircraft in ('myrocket','myrocket2','myrocket_cluster','myrocket_series'):
        path=root/'aircraft'/aircraft/(aircraft+'.xml');s=path.read_text(encoding='utf-8')
        def function(match):
            block=match.group();name=re.search(r'<function name="([^"]+)"',block).group(1).split('/')[-1]
            kind=next(k for k in KINDS if name==k or name==k+'_sust')
            c=COEFF_STACK if aircraft=='myrocket2' and not name.endswith('_sust') else COEFF_SINGLE
            def table(m):
                angle=re.search(r'<independentVar[^>]*>(.*?)</independentVar>',m.group(),re.S).group(1).strip()
                if angle not in ('aero/alpha-rad','aero/beta-rad'): return m.group()
                if kind in ('Lift_alpha','Side_beta'):
                    sign=.5 if kind=='Lift_alpha' else -.5
                    return f'<product><value>{sign}</value><sin><product><value>2</value><property>{angle}</property></product></sin></product>'
                if kind in ('Drag_crossflow','Drag_beta'):
                    return f'<product><value>{c:.15g}</value>'+''.join(f'<abs><sin><property>{angle}</property></sin></abs>' for _ in range(3))+'</product>'
                sign=1 if kind=='Lift_crossflow' else -1
                return f'<product><value>{sign*c:.15g}</value><sin><property>{angle}</property></sin><abs><sin><property>{angle}</property></sin></abs><cos><property>{angle}</property></cos></product>'
            return re.sub(r'<table>.*?</table>',table,block,flags=re.S)
        s=re.sub(r'<function name="aero/force/(?:'+ '|'.join(KINDS)+r')(?:_sust)?">.*?</function>',function,s,flags=re.S)
        path.write_text(s,encoding='utf-8')


def verify(root):
    records=[]
    for aircraft in ('myrocket','myrocket2','myrocket_cluster','myrocket_series'):
        path=root/'aircraft'/aircraft/(aircraft+'.xml');xml=ET.parse(path).getroot()
        funcs={f.attrib['name'].split('/')[-1]:f for f in xml.findall('aerodynamics/axis/function')}
        phases=('', '_sust') if aircraft=='myrocket2' else ('',)
        for suffix in phases:
            for kind in KINDS:
                node=funcs[kind+suffix]
                assert not any(t.findtext('independentVar') in ('aero/alpha-rad','aero/beta-rad') for t in node.findall('.//table')), (aircraft,kind)
        # Both signs and finite angles, including the small-alpha derivative.
        def evaluate(node, props):
            tag=node.tag;children=list(node)
            if tag=='value':return float(node.text)
            if tag=='property':return props[node.text.strip()]
            if tag in ('function','product'):return math.prod(evaluate(c,props) for c in children if c.tag!='description')
            if tag=='sin':return math.sin(evaluate(children[0],props))
            if tag=='cos':return math.cos(evaluate(children[0],props))
            if tag=='abs':return abs(evaluate(children[0],props))
            if tag=='sum':return sum(evaluate(c,props) for c in children)
            if tag=='difference':return evaluate(children[0],props)-evaluate(children[1],props)
            if tag=='table':
                x=props[node.findtext('independentVar').strip()]
                points=[tuple(map(float,l.split())) for l in node.findtext('tableData').splitlines() if l.strip()]
                if x<=points[0][0]:return points[0][1]
                for (x0,y0),(x1,y1) in zip(points,points[1:]):
                    if x<=x1:return y0+(y1-y0)*(x-x0)/(x1-x0)
                return points[-1][1]
            raise ValueError(tag)
        maxerr=0.
        for suffix in phases:
            c=COEFF_STACK if aircraft=='myrocket2' and suffix=='' else COEFF_SINGLE
            stage=1 if suffix else 0
            for a in (-.26,-.0174532925199433,-1e-7,0.,1e-7,.0174532925199433,.26):
                props={'fcs/body-aero-gain':1.,'fcs/k-cn':1.,'fcs/k-cd':1.,'fcs/off-rail':1.,'fcs/staged':stage,
                       'aero/qbar-psf':1.,'metrics/Sw-sqft':1.,'velocities/mach':0.,'aero/alpha-rad':a,'aero/beta-rad':a}
                for kind in KINDS:
                    node=funcs[kind+suffix];got=evaluate(node,props)
                    if kind in ('Lift_alpha','Side_beta'):
                        cn=evaluate(node.find('.//table'),props);expected=(1 if kind=='Lift_alpha' else -1)*cn*.5*math.sin(2*a)
                    elif kind in ('Drag_crossflow','Drag_beta'):expected=c*abs(math.sin(a))**3
                    else:expected=(1 if kind=='Lift_crossflow' else -1)*c*math.sin(a)*abs(math.sin(a))*math.cos(a)
                    maxerr=max(maxerr,abs(got-expected));assert math.isclose(got,expected,rel_tol=1e-12,abs_tol=1e-12),(aircraft,kind,a,got,expected)
        f=jsbsim.FGFDMExec(str(root));f.set_debug_level(0)
        assert f.load_model(aircraft);assert f.load_ic('reset',True);assert f.run_ic();f.disable_output();assert f.run()
        records.append(dict(aircraft=aircraft,phases=len(phases),forces_checked=len(phases)*6,
                            angle_tables_remaining=0,max_expression_error=maxerr,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    assert sum(r['forces_checked'] for r in records)==30
    return records


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',default=str(Path(__file__).resolve().parent));parser.add_argument('--patch',action='store_true');parser.add_argument('--json',required=True);args=parser.parse_args()
    root=Path(args.root)
    if args.patch:patch(root)
    records=verify(root);Path(args.json).write_text(json.dumps(records,indent=2));print(json.dumps(records,indent=2))
