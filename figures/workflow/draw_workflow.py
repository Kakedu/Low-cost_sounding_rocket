from pathlib import Path
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description="Render the model and analysis workflow diagram.")
    parser.add_argument("--output-dir", type=Path, default=HERE / "generated",
                        help="Figure destination (default: generated beside this script)")
    args=parser.parse_args()
    output=args.output_dir.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    fig,ax=plt.subplots(figsize=(11.2,8.7));fig.subplots_adjust(0,0,1,1);ax.set(xlim=(0,11.2),ylim=(0,8.7));ax.axis('off')
    edge='#adb3b7';blue='#2b6cb0';ink='#202529'
    def box(x,y,w,h,title,lines,hi=False,fs=13.7):
     p=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.07,rounding_size=0.10',linewidth=1.25,edgecolor=blue if hi else edge,facecolor='#edf4fa' if hi else 'white');ax.add_patch(p)
     ax.text(x+.11,y+h-.16,title,fontsize=15.0,fontweight='bold',va='top',color=ink)
     ax.text(x+.11,y+h-.47,'\n'.join(lines),fontsize=fs,va='top',linespacing=1.3,color='#4b5155')
    def arrow(x1,y1,x2,y2):ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=13,linewidth=1.25,color='#929ca3'))
    ax.text(.25,8.42,'Model and analysis workflow',fontsize=21,fontweight='bold',color=ink)
    ax.text(.35,7.92,'Vehicle definition in XML',fontsize=16,fontweight='bold',color=ink)
    box(.30,6.90,3.45,.80,'Mass and geometry',['mass, CG, inertia, CP'])
    box(.30,5.91,3.45,.80,'Propulsion',['thrust versus burned mass'])
    box(.30,4.92,3.45,.80,'Aerodynamics',['drag, normal force, crossflow'])
    box(.30,3.93,3.45,.80,'Flight-control channels',['rail, brake, recovery, staging'],True)
    box(.30,2.94,3.45,.80,'Launch and initial conditions',['ignition, attitude, flight events'])
    box(4.05,4.80,1.90,1.05,'JSBSim 1.3.1',['6-DOF motion','1 ms time step'])
    arrow(3.82,5.26,3.98,5.26)
    ax.text(7.03,7.92,'Python study drivers',fontsize=16,fontweight='bold',color=ink)
    box(7.00,6.57,3.84,1.10,'sweep.py',['one-at-a-time, Morris, Monte Carlo,','brake range and design map'])
    box(7.00,5.13,3.84,1.17,'sil.py',['sensors and state estimation,','apogee prediction, brake commands','50 Hz flight computer'],True)
    box(7.00,3.94,3.84,.90,'anglesweep.py',['launch-angle and wind comparisons'])
    box(7.00,2.50,3.84,1.14,'Staging studies',['stage_compare.py, sweep_staging.py,','staging_seq.py, staging_post.py'])
    arrow(6.93,5.65,6.02,5.65);arrow(6.02,4.99,6.93,4.99)
    ax.text(6.48,6.05,'inputs /\ncommands',fontsize=10.5,ha='center',va='bottom',color='#4b5155')
    ax.text(6.48,4.75,'states /\nresults',fontsize=10.5,ha='center',va='top',color='#4b5155')
    box(.30,.34,5.05,1.42,'Numerical checks  Section 3.10',['crosscheck.py, crosscheck_staging.py','Independent point-mass equations','Same vehicle and environment inputs'])
    box(5.72,.34,5.12,1.42,'Flight comparisons  Section 3.11',['Prometheus: GPS and barometric records','Cavour: barometric record','Separate models using published vehicle data'])
    arrow(4.67,4.72,4.67,1.83);arrow(5.55,4.72,6.10,1.83)
    fig.savefig(output/'model_workflow.png',dpi=210,facecolor='white');fig.savefig(output/'model_workflow.svg',facecolor='white');plt.close(fig)
    print('Workflow diagram saved')


if __name__ == "__main__":
    main()
