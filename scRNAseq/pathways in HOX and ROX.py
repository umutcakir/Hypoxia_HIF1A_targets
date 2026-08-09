"""
Figure 1F  -  pathway-level shift of glycolysis and oxidative phosphorylation
              under hypoxia and roxadustat, in monocytes and CD4+ T cells.

INPUT   S7_all_celltypes_all_comparisons_bimod_all_results.csv
        (single-cell differential expression: gene, cell_type, comparison,
         avg_log2FC, p_val_adj)
OUTPUT  Figure1F_pathway_shift_Mono_CD4T_vertical.png  (1200 dpi)
        Figure1F_pathway_shift_Mono_CD4T_vertical.svg  (vector)

METHOD  Gene sets are defined by pathway membership among all genes TESTED in
        that cell type (no significance pre-selection, which would bias the
        hypoxia vs roxadustat comparison). Pseudogenes, antisense and
        readthrough transcripts are excluded by symbol pattern; mitochondrial
        genes are exempt from that filter so MT-ATP6 / MT-ATP8 are retained.
        Each violin is tested against the distribution of ALL tested genes in
        the same cell type and condition (two-sided Mann-Whitney), so the grey
        band (interquartile range of all tested genes) is the reference, not 0.
        Panels use independent y-scales.
"""
DATA = "S7_all_celltypes_all_comparisons_bimod_all_results.csv"   # <-- set path

import pandas as pd, numpy as np, re, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.text as mtext
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu
import warnings; warnings.filterwarnings('ignore')
BADRE=re.compile(r'(P\d+$|-AS\d*$|-DT$|-IT\d*$|-KCTD14$)')
def _ok(g): return g.startswith('MT-') or not BADRE.search(g)   # never filter mtDNA genes
plt.rcParams.update({'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black'})

d=pd.read_csv(DATA)
glyc=['HK1','HK2','HK3','GPI','PFKL','PFKM','PFKP','ALDOA','ALDOB','ALDOC','TPI1','GAPDH','PGK1',
      'PGAM1','ENO1','ENO2','ENO3','PKM','LDHA','SLC2A1','SLC2A3','SLC16A3']
def members(tested):
    return [('Glycolysis',[g for g in glyc if g in tested],'#c0392b'),
            ('Complex I',[g for g in tested if re.match(r'NDUF',g) and _ok(g)],'#2c6fa8'),
            ('Complex II/III',[g for g in tested if re.match(r'(SDH|UQCR)',g) and _ok(g)],'#2c6fa8'),
            ('Complex IV',[g for g in tested if re.match(r'COX',g) and not re.match(r'COX1[0-9]|COX4I2',g) and _ok(g)],'#2c6fa8'),
            ('ATP synthase',[g for g in tested if re.match(r'ATP5',g) and _ok(g)],'#2c6fa8'),
            ('OXPHOS\n(mtDNA)',[g for g in tested if re.match(r'MT-(CO|ND|ATP|CYB)',g) and _ok(g)],'#1a4f7a')]
conds=[('HOX_vs_NOX','Hyp',-0.20,0.62),('ROX_vs_NOX','Rox',0.20,0.28)]
panels=[('Monocytes','a'),('CD4+ T cells','b')]

fig,axes=plt.subplots(2,1,figsize=(8.4,11.2))
fig.subplots_adjust(left=0.115,right=0.99,top=0.945,bottom=0.085,hspace=0.46)
ylims={'Monocytes':(-1.25,3.35),'CD4+ T cells':(-0.75,2.15)}

for ax,(ct,tag) in zip(axes,panels):
    mono={c:d[(d.comparison==c)&(d.cell_type==ct)] for c,_,_,_ in conds}
    bgH=mono['HOX_vs_NOX'].avg_log2FC.dropna().values
    ax.axhspan(np.percentile(bgH,25),np.percentile(bgH,75),color='#dddddd',alpha=0.6,zorder=0)
    ax.axhline(0,color='#555555',lw=1.1,zorder=1)
    tr=ax.get_xaxis_transform()
    labels=[s[0] for s in members(set(mono['HOX_vs_NOX'].gene))]
    for comp,cname,off,alpha in conds:
        m=mono[comp]; bg=m.avg_log2FC.dropna().values
        for i,(lab,gl,col) in enumerate(members(set(m.gene))):
            v=m[m.gene.isin(gl)].avg_log2FC.dropna().values
            if len(v)==0: continue
            p=mannwhitneyu(v,bg,alternative='two-sided').pvalue
            parts=ax.violinplot([v],positions=[i+off],widths=0.34,showextrema=False)
            for pc in parts['bodies']:
                pc.set_facecolor(col); pc.set_alpha(alpha); pc.set_edgecolor('black'); pc.set_linewidth(0.6)
            ax.scatter(np.random.normal(i+off,0.035,len(v)),v,s=7,color=col,edgecolors='black',
                       linewidths=0.22,alpha=0.85,zorder=3)
            ax.plot([i+off-0.15,i+off+0.15],[np.median(v)]*2,color='black',lw=1.8,zorder=4)
            ax.text(i,1.115 if comp=='HOX_vs_NOX' else 1.055,f'{cname} P = {p:.1e}',transform=tr,
                    ha='center',va='bottom',fontsize=7.2,clip_on=False)
            if comp=='HOX_vs_NOX':
                ax.text(i,1.012,f'n = {len(v)}',transform=tr,ha='center',va='bottom',fontsize=7.2,clip_on=False)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels,fontsize=8.4,rotation=25,ha='right')
    ax.set_ylim(*ylims[ct]); ax.set_xlim(-0.6,len(labels)-0.4)
    ax.grid(True,axis='y',color='#eeeeee',lw=0.6); ax.set_axisbelow(True)
    for sp in ax.spines.values(): sp.set_color('#999999'); sp.set_linewidth(0.7)
    ax.text(len(labels)-0.45,np.percentile(bgH,75)+0.05,'IQR of all tested genes',
            fontsize=6.8,ha='right',va='bottom',color='black')

for a_,(ct_,_) in zip(axes,panels):
    a_.set_ylabel('log$_2$ fold change vs normoxia\n(%s)'%ct_.replace('CD4+','CD4$^+$'),fontsize=9.8)
axes[0].legend(handles=[Patch(facecolor='#777777',alpha=0.62,edgecolor='black',label='Hypoxia'),
                        Patch(facecolor='#777777',alpha=0.28,edgecolor='black',label='Roxadustat')],
               loc='upper right',frameon=False,fontsize=8.2,handlelength=1.5,bbox_to_anchor=(1.0,0.985))
for o in fig.findobj(mtext.Text): o.set_color('black')
fig.savefig('Figure1F_pathway_shift_Mono_CD4T_vertical.png',dpi=1200,facecolor='white')
fig.savefig('Figure1F_pathway_shift_Mono_CD4T_vertical.svg',facecolor='white')





"""
Supplementary figure  -  glycolysis and OXPHOS across all PBMC cell types.

  a  pathway level: one point per gene set per cell type.
     colour = median log2 fold change of all tested genes in the set
     size   = -log10 P, two-sided Mann-Whitney vs all tested genes in that cell type
  b  individual genes of the same six sets, split into two column blocks.
     colour = log2 fold change, size = -log10 p_adj of the per-gene test,
     open circles = not significant, pale grey = not tested in that cell type.

INPUT   S7_all_celltypes_all_comparisons_bimod_all_results.csv
OUTPUT  SupplementaryFigure_OXPHOS_combined.png  (600 dpi)
        SupplementaryFigure_OXPHOS_combined.svg  (vector)

NOTE    Gene sets are identical to those used in Figure 1F. pDC is absent
        because the dataset contains no HOX_vs_NOX rows for that population.
"""
DATA = "S7_all_celltypes_all_comparisons_bimod_all_results.csv"   # <-- set path

import pandas as pd, numpy as np, re, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.text as mtext
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from matplotlib.transforms import blended_transform_factory
from scipy.stats import mannwhitneyu
import warnings; warnings.filterwarnings('ignore')
BADRE=re.compile(r'(P\d+$|-AS\d*$|-DT$|-IT\d*$|-KCTD14$)')
def _ok(g): return g.startswith('MT-') or not BADRE.search(g)
plt.rcParams.update({'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black'})

d=pd.read_csv(DATA)
cts=['Monocytes','CD4+ T cells','NK cells','B cells','MAIT cells','Neutrophils','Basophils','Proliferating cells']
CTL=[c.replace('Proliferating cells','Proliferating').replace('CD4+','CD4$^+$') for c in cts]
comps=[('HOX_vs_NOX','Hypoxia'),('ROX_vs_NOX','Roxadustat')]
glyc=['HK1','HK2','HK3','GPI','PFKL','PFKM','PFKP','ALDOA','ALDOB','ALDOC','TPI1','GAPDH','PGK1',
      'PGAM1','ENO1','ENO2','ENO3','PKM','LDHA','SLC2A1','SLC2A3','SLC16A3']
def members(t):
    return [('Glycolysis',[g for g in glyc if g in t]),
            ('Complex I',[g for g in t if re.match(r'NDUF',g) and _ok(g)]),
            ('Complex II/III',[g for g in t if re.match(r'(SDH|UQCR)',g) and _ok(g)]),
            ('Complex IV',[g for g in t if re.match(r'COX',g) and not re.match(r'COX1[0-9]|COX4I2',g) and _ok(g)]),
            ('ATP synthase',[g for g in t if re.match(r'ATP5',g) and _ok(g)]),
            ('OXPHOS (mtDNA)',[g for g in t if re.match(r'MT-(CO|ND|ATP|CYB)',g) and _ok(g)])]
setnames=[s[0] for s in members(set(d.gene))]
T=set(d[d.cell_type.isin(cts)].gene)
LABMAP={'Glycolysis':'Glycolysis','Complex I':'Complex I\n(nuclear)','Complex II/III':'Complex II/III\n(nuclear)',
        'Complex IV':'Complex IV\n(nuclear)','ATP synthase':'ATP synthase\n(nuclear)','OXPHOS (mtDNA)':'OXPHOS\n(mtDNA)'}
G=[(LABMAP[l],gl) for l,gl in members(T)]
h=d[(d.comparison=='HOX_vs_NOX')&d.cell_type.isin(cts)].copy()
h['nlp']=-np.log10(h.p_val_adj.clip(lower=1e-100)); rank=h.groupby('gene')['nlp'].mean()
BG=[]
for bl in [G[:2],G[2:]]:
    genes=[];bk=[]
    for lab,gl in bl:
        gl=sorted(gl,key=lambda g:-rank.get(g,0)); bk.append((lab,len(genes),len(genes)+len(gl))); genes+=gl
    BG.append((genes,bk))
nmax=max(len(g) for g,_ in BG); vmax=1.0
cmap=plt.get_cmap('RdBu_r'); norm=plt.Normalize(-vmax,vmax)
sz  =lambda p: float(np.clip(-np.log10(max(p,1e-100)),0,40))*1.8+9
szS =lambda p: float(np.clip(-np.log10(max(p,1e-30)),0,20))*13+22

H=14.1
fig=plt.figure(figsize=(15.2,H))
# ---------- panel a : set level ----------
ya,ha=0.895,0.082
for k,(comp,title) in enumerate(comps):
    ax=fig.add_axes([0.105+k*0.198,ya,0.176,ha])
    for j,ct in enumerate(cts):
        m=d[(d.comparison==comp)&(d.cell_type==ct)]
        if len(m)==0: continue
        bg=m.avg_log2FC.dropna().values
        for i,(lab,gl) in enumerate(members(set(m.gene))):
            v=m[m.gene.isin(gl)].avg_log2FC.dropna().values
            if len(v)==0: continue
            p=mannwhitneyu(v,bg,alternative='two-sided').pvalue
            y=len(setnames)-1-i
            if p<0.05: ax.scatter(j,y,s=szS(p),c=[cmap(norm(np.median(v)))],edgecolors='black',lw=0.5,zorder=3)
            else: ax.scatter(j,y,s=30,facecolors='none',edgecolors='#999999',lw=0.8,zorder=2)
    ax.set_xticks(range(len(cts))); ax.set_xticklabels(CTL,rotation=90,ha='center',va='top',fontsize=7.6)
    ax.set_xlim(-0.6,len(cts)-0.4); ax.set_ylim(-0.7,len(setnames)-0.3)
    ax.set_title(title,fontsize=9.6,pad=6)
    ax.grid(True,color='#eeeeee',lw=0.6,zorder=0); ax.set_axisbelow(True)
    for sp in ax.spines.values(): sp.set_color('#999999'); sp.set_linewidth(0.7)
    if k==0:
        ax.set_yticks(range(len(setnames))); ax.set_yticklabels(setnames[::-1],fontsize=8.2)
        ax.text(-0.52,1.09,'a',transform=ax.transAxes,fontsize=13,fontweight='bold',va='bottom')
    else: ax.set_yticks([])
fig.text(0.105,0.038,'a, pathway level: colour, median log$_2$ fold change of all tested genes in the set; size, $-$log$_{10}$ P from a two-sided',fontsize=8.0,ha='left')
fig.text(0.105,0.019,'Mann-Whitney test against all tested genes in that cell type.  b, individual genes of the same six sets.',fontsize=8.0,ha='left')
# ---------- panel b : individual genes ----------
gs=GridSpec(1,5,width_ratios=[1,1,0.42,1,1],wspace=0.07,left=0.105,right=0.865,top=0.800,bottom=0.135)
for bi,(genes,bk) in enumerate(BG):
    n=len(genes)
    for ci,(comp,title) in enumerate(comps):
        ax=fig.add_subplot(gs[0,(0 if bi==0 else 3)+ci])
        sub=d[(d.comparison==comp)&d.cell_type.isin(cts)]
        for j,ct in enumerate(cts):
            s=sub[sub.cell_type==ct].set_index('gene')
            for i,g in enumerate(genes):
                y=nmax-1-i
                if g not in s.index:
                    ax.scatter(j,y,s=7,facecolors='none',edgecolors='#cccccc',lw=0.35,zorder=2); continue
                r=s.loc[g]
                if r.p_val_adj<0.05:
                    ax.scatter(j,y,s=sz(r.p_val_adj),c=[cmap(norm(r.avg_log2FC))],edgecolors='black',lw=0.3,zorder=3)
                else:
                    ax.scatter(j,y,s=13,facecolors='none',edgecolors='#b0b0b0',lw=0.5,zorder=2)
        ax.set_xticks(range(len(cts))); ax.set_xticklabels(CTL,rotation=90,ha='center',va='top',fontsize=7.2)
        ax.set_xlim(-0.6,len(cts)-0.4); ax.set_ylim(-0.7,nmax-0.3)
        ax.set_title(title,fontsize=9.4,pad=6)
        ax.grid(True,color='#f2f2f2',lw=0.45,zorder=0); ax.set_axisbelow(True)
        for sp in ax.spines.values(): sp.set_color('#999999'); sp.set_linewidth(0.6)
        for _,s0,e0 in bk[:-1]: ax.axhline(nmax-1-e0+0.5,color='#666666',lw=0.7,ls=(0,(4,3)),zorder=1)
        if ci==0:
            ax.set_yticks(range(nmax-n,nmax)); ax.set_yticklabels(genes[::-1],fontsize=5.6)
            bt=blended_transform_factory(ax.transAxes,ax.transData)
            for lab,s0,e0 in bk:
                yt=nmax-1-s0+0.45; yb=nmax-1-(e0-1)-0.45
                ax.plot([-0.245,-0.245],[yb,yt],transform=bt,color='black',lw=1.1,clip_on=False)
                ax.text(-0.315,(yt+yb)/2,lab,transform=bt,rotation=90,va='center',ha='center',
                        fontsize=7.4,color='black',clip_on=False)
            if bi==0: ax.text(-0.40,1.035,'b',transform=ax.transAxes,fontsize=13,fontweight='bold',va='bottom')
        else: ax.set_yticks([])
# ---------- shared legends ----------
cax=fig.add_axes([0.878,0.60,0.011,0.15])
cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
cb.set_ticks([-1,0,1]); cb.set_label('log$_2$ fold change',fontsize=8.4); cb.ax.tick_params(labelsize=7.5)
hs=[Line2D([],[],marker='o',color='none',markerfacecolor='#888888',markeredgecolor='black',
    markersize=np.sqrt(sz(10**-v)),label=f'{v}') for v in [2,10,20,40]]
hs.append(Line2D([],[],marker='o',color='none',markerfacecolor='none',markeredgecolor='#b0b0b0',
    markersize=np.sqrt(13),label='n.s.'))
l1=fig.legend(handles=hs,title='$-$log$_{10}$ p$_{adj}$ (gene)',loc='upper left',bbox_to_anchor=(0.870,0.55),
    frameon=False,fontsize=7.6,title_fontsize=8.2,labelspacing=0.95,handletextpad=1.1)
hs2=[Line2D([],[],marker='o',color='none',markerfacecolor='#888888',markeredgecolor='black',
     markersize=np.sqrt(szS(10**-v)),label=f'{v}') for v in [2,5,10,20]]
hs2.append(Line2D([],[],marker='o',color='none',markerfacecolor='none',markeredgecolor='#999999',
     markersize=np.sqrt(30),label='n.s.'))
l2=fig.legend(handles=hs2,title='$-$log$_{10}$ P (pathway)',loc='upper left',bbox_to_anchor=(0.870,0.86),
     frameon=False,fontsize=7.6,title_fontsize=8.2,labelspacing=1.0,handletextpad=1.1)
for o in fig.findobj(mtext.Text): o.set_color('black')
l1.get_title().set_color('black'); l2.get_title().set_color('black')
fig.savefig('SupplementaryFigure_OXPHOS_combined.png',dpi=600,facecolor='white')
fig.savefig('SupplementaryFigure_OXPHOS_combined.svg',facecolor='white')
print("ok  blocks:",[len(g) for g,_ in BG])
print("saved")
