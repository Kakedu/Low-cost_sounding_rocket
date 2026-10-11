# Converts the RASP .eng thrust curve into a JSBSim engine table.
# JSBSim looks thrust up by propellant burned, not by time.

LBF=4.4482216; LB=0.45359237
eng=[l.split() for l in open("M1670.eng") if l.strip() and not l.startswith(";")]
hdr=eng[0]; pts=[(0.0,0.0)]+[(float(a),float(b)) for a,b in eng[1:]]
prop_kg=float(hdr[4]); total_kg=float(hdr[5])
prop_lb=prop_kg/LB
I=0.0; cum=[0.0]
for (t0,f0),(t1,f1) in zip(pts,pts[1:]):
    I+=0.5*(f0+f1)*(t1-t0); cum.append(I)
I_lbfs=I/LBF; isp=I_lbfs/prop_lb
print(f"impulse {I:.1f} N.s  prop {prop_kg*1000:.0f} g = {prop_lb:.4f} lb  Isp {isp:.2f} s  burnout mass {(total_kg-prop_kg)*1000:.0f} g")
rows=[(prop_lb*c/I, f/LBF, t) for (t,f),c in zip(pts,cum)]
first=rows[1]; build=first[2]
table=[(0.0,first[1],first[2])]+rows[2:-1]
last_nonzero=rows[-2][1]
table.append((prop_lb, rows[-1][1] if rows[-1][1]>0 else last_nonzero*0.5, rows[-1][2]))
lines="\n".join(f"      {p:9.4f}  {F:9.2f}   <!-- t = {t:5.3f} s -->" for p,F,t in table)
xml=f'''<?xml version="1.0"?>
<!-- Cesaroni 6026M1670-P, {I:.0f} N.s, Isp {isp:.1f} s. Thrust indexed by propellant burned (lb). -->
<rocket_engine name="M1670">

  <isp> {isp:.3f} </isp>
  <builduptime> {build:.3f} </builduptime>

  <thrust_table name="propulsion/thrust_prop_remain" type="internal">
    <tableData>
{lines}
    </tableData>
  </thrust_table>

</rocket_engine>
'''
open("M1670_new.xml","w").write(xml)
