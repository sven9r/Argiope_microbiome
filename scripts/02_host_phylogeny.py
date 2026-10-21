#!/usr/bin/env python3
"""Render species read-composition pies and detection pies on published host tree."""
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Patch, Wedge, Polygon, Rectangle
from matplotlib.colors import to_rgb
from Bio import Phylo


# Simplified national flags for the recorded specimen origins.
LOCALITIES = {
    'Australia': 'AU', 'Brazil': 'BR', 'Cape Verde': 'CV', 'El Salvador': 'SV',
    'France': 'FR', 'French Guiana': 'GF', 'Germany': 'DE', 'Greece': 'GR',
    'India': 'IN', 'Indonesia': 'ID', 'Ivory Coast': 'CI', 'Japan': 'JP',
    'Laos': 'LA', 'Madagascar': 'MG', 'Malaysia': 'MY', 'Micronesia': 'FM',
    'Myanmar': 'MM', 'Papua New Guinea': 'PG', 'Philippines': 'PH',
    'Russia': 'RU', 'Senegal': 'SN', 'Spain': 'ES', 'Tanzania': 'TZ', 'USA': 'US',
}


def draw_flag(ax, country, x, y, width=.80, height=.48, eligible=True):
    """Draw a centred icon, locality code, and dagger for no eligible profile."""
    code = LOCALITIES[country]
    left, top = x-width/2, y-height/2-.12
    border = Rectangle((left, top), width, height, facecolor='white',
                       edgecolor='none', zorder=3)
    ax.add_patch(border)

    def rect(u, v, w, h, color):
        patch=Rectangle((left+u*width, top+v*height), w*width, h*height,
                        facecolor=color, edgecolor='none', zorder=4)
        patch.set_clip_path(border);ax.add_patch(patch)

    def poly(points, color):
        patch=Polygon([(left+u*width, top+v*height) for u,v in points],
                      closed=True, facecolor=color, edgecolor='none', zorder=5)
        patch.set_clip_path(border);ax.add_patch(patch)

    def circle(u, v, radius, color):
        patch=Circle((left+u*width, top+v*height), radius*height,
                     facecolor=color, edgecolor='none', zorder=6)
        patch.set_clip_path(border);ax.add_patch(patch)

    def star(u, v, radius, color, points=5):
        # Correct x for the rectangular flag so stars remain circular.
        xy=[]
        for n in range(points*2):
            a=-math.pi/2+n*math.pi/points
            r=radius if n%2==0 else radius*.42
            xy.append((u+math.cos(a)*r*height/width, v+math.sin(a)*r))
        poly(xy,color)

    def horizontal(colors, sizes=None):
        sizes=sizes or [1/len(colors)]*len(colors)
        v=0
        for color,size in zip(colors,sizes):rect(0,v,1,size,color);v+=size

    def vertical(colors):
        for n,color in enumerate(colors):rect(n/len(colors),0,1/len(colors),1,color)

    def union_jack(w=.48,h=.48):
        rect(0,0,w,h,'#012169')
        poly([(0,0),(w*.13,0),(w,h*.80),(w,h),(w*.87,h),(0,h*.20)],'white')
        poly([(0,h),(0,h*.80),(w*.87,0),(w,0),(w,h*.20),(w*.13,h)],'white')
        poly([(0,0),(w*.055,0),(w,h*.91),(w,h),(w*.945,h),(0,h*.09)],'#c8102e')
        poly([(0,h),(0,h*.91),(w*.945,0),(w,0),(w,h*.09),(w*.055,h)],'#c8102e')
        rect(w*.40,0,w*.20,h,'white');rect(0,h*.36,w,h*.28,'white')
        rect(w*.445,0,w*.11,h,'#c8102e');rect(0,h*.425,w,h*.15,'#c8102e')

    if code in ('FR','GF'):vertical(['#002395','white','#ed2939'])
    elif code=='DE':horizontal(['black','#dd0000','#ffce00'])
    elif code=='RU':horizontal(['white','#0039a6','#d52b1e'])
    elif code=='ID':horizontal(['#ff0000','white'])
    elif code=='JP':circle(.5,.5,.30,'#bc002d')
    elif code=='CI':vertical(['#f77f00','white','#009e60'])
    elif code=='SN':
        vertical(['#00853f','#fdef42','#e31b23']);star(.5,.5,.19,'#00853f')
    elif code=='MG':
        rect(1/3,0,2/3,.5,'#fc3d32');rect(1/3,.5,2/3,.5,'#007e3a')
    elif code=='LA':
        horizontal(['#ce1126','#002868','#ce1126'],[.25,.5,.25]);circle(.5,.5,.20,'white')
    elif code=='MM':
        horizontal(['#fecb00','#34b233','#ea2839']);star(.5,.5,.32,'white')
    elif code=='IN':
        horizontal(['#ff9933','white','#138808']);circle(.5,.5,.135,'#000080');circle(.5,.5,.10,'white')
        for a in range(24):
            rad=2*math.pi*a/24
            ax.plot([left+.5*width,left+.5*width+.13*height*math.cos(rad)],
                    [top+.5*height,top+.5*height+.13*height*math.sin(rad)],
                    color='#000080',lw=.20,zorder=7)
    elif code=='BR':
        rect(0,0,1,1,'#009739');poly([(.08,.5),(.5,.08),(.92,.5),(.5,.92)],'#ffdf00')
        circle(.5,.5,.27,'#012169');poly([(.35,.43),(.66,.54),(.65,.60),(.35,.49)],'white')
    elif code=='SV':
        horizontal(['#0047ab','white','#0047ab']);poly([(.44,.61),(.5,.39),(.56,.61)],'#b9a861')
    elif code=='ES':
        horizontal(['#aa151b','#f1bf00','#aa151b'],[.25,.5,.25])
        rect(.25,.38,.11,.24,'#aa151b');rect(.25,.38,.055,.11,'#f1bf00');rect(.31,.49,.05,.12,'white')
    elif code=='GR':
        horizontal(['#0d5eaf' if n%2==0 else 'white' for n in range(9)])
        rect(0,0,.38,5/9,'#0d5eaf');rect(.15,0,.075,5/9,'white');rect(0,2/9,.38,1/9,'white')
    elif code=='CV':
        rect(0,0,1,1,'#003893');rect(0,.50,1,.10,'white');rect(0,.60,1,.10,'#cf2027');rect(0,.70,1,.10,'white')
        for n in range(10):
            a=2*math.pi*n/10-math.pi/2
            star(.33+math.cos(a)*.17,.65+math.sin(a)*.28,.035,'#f7d116')
    elif code=='TZ':
        rect(0,0,1,1,'#1eb53a');poly([(0,1),(1,0),(1,1)],'#00a3dd')
        poly([(0,.77),(.77,0),(1,0),(1,.23),(.23,1),(0,1)],'#fcd116')
        poly([(0,.84),(.84,0),(1,0),(1,.16),(.16,1),(0,1)],'black')
    elif code=='US':
        horizontal(['#b22234' if n%2==0 else 'white' for n in range(13)])
        rect(0,0,.40,7/13,'#3c3b6e')
        for r in range(9):
            for c in range(6 if r%2==0 else 5):
                circle(.032+c*.064+(.032 if r%2 else 0),.028+r*.060,.010,'white')
    elif code=='AU':
        rect(0,0,1,1,'#012169');union_jack()
        star(.24,.75,.13,'white',7)
        for u,v,r,p in [(.72,.23,.09,7),(.61,.54,.08,7),(.85,.47,.08,7),(.75,.82,.09,7),(.77,.60,.045,5)]:
            star(u,v,r,'white',p)
    elif code=='MY':
        horizontal(['#cc0000' if n%2==0 else 'white' for n in range(14)])
        rect(0,0,.50,8/14,'#010066');circle(.22,.28,.205,'#ffcc00');circle(.25,.255,.17,'#010066');star(.36,.27,.16,'#ffcc00',14)
    elif code=='PH':
        horizontal(['#0038a8','#ce1126']);poly([(0,0),(.48,.5),(0,1)],'white')
        circle(.15,.5,.095,'#fcd116')
        for a in range(8):
            rad=2*math.pi*a/8
            star(.15+math.cos(rad)*.11*height/width,.5+math.sin(rad)*.11,.032,'#fcd116',4)
        for u,v in [(.035,.11),(.035,.89),(.39,.5)]:star(u,v,.055,'#fcd116')
    elif code=='PG':
        rect(0,0,1,1,'black');poly([(0,0),(1,0),(1,1)],'#ce1126')
        # Bird-of-paradise silhouette and the Southern Cross.
        poly([(.62,.18),(.69,.33),(.89,.16),(.78,.40),(.86,.61),(.70,.43),(.57,.44),(.68,.35),(.55,.28)],'#fcd116')
        for u,v,r in [(.23,.40,.055),(.13,.62,.045),(.34,.59,.045),(.24,.80,.055),(.30,.71,.028)]:star(u,v,r,'white')
    elif code=='FM':
        rect(0,0,1,1,'#75b2dd')
        for u,v in [(.5,.20),(.32,.5),(.68,.5),(.5,.80)]:star(u,v,.105,'white')
    else:raise ValueError(f'No recorded-country flag definition: {country}')

    ax.add_patch(Rectangle((left,top),width,height,fill=False,edgecolor='#78848b',
                           linewidth=.50,linestyle='--' if not eligible else '-',zorder=8))
    ax.text(x,y+.28,code+('†' if not eligible else ''),ha='center',va='center',
            fontsize=6.8,color='#243641',zorder=9)


ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output/species_symbionts'
COLORS={'DUSA-like OTU2':'#ff6a0a','DUSA close sequence relatives':'#efb17e','Mesoplasma':'#eeee00',
        'Wolbachia':'#9067b6','Spiroplasma':'#318c91','Other named candidates':'#c04b5a',
        'Unresolved candidate annotations':'#6f91b0','Other 16S features':'#d8d8d8'}
KEYS=list(COLORS)
NAMES={'DUSA-like OTU2':'DUSA-like OTU2','DUSA close sequence relatives':'Close to DUSA (OTU287)',
       'Mesoplasma':'Mesoplasma','Wolbachia':'Wolbachia','Spiroplasma':'Spiroplasma',
       'Other named candidates':'Other named candidate genera','Unresolved candidate annotations':'Unresolved candidate annotations',
       'Other 16S features':'Other 16S features'}

def read(path):
    with path.open() as f:return list(csv.DictReader(f))

def bucket(group):
    if group in KEYS:return group
    if any(x in group.lower() for x in ('hint','unresolved','unnamed','nonbacterial')):
        return 'Unresolved candidate annotations'
    return 'Other named candidates'

def pie(ax,x,y,values,colors,r=.31):
    start=90
    for value,color in zip(values,colors):
        if value<=0:continue
        end=start+360*value
        ax.add_patch(Wedge((x,y),r,start,end,facecolor=color,edgecolor='white',linewidth=.35))
        start=end
    ax.add_patch(Circle((x,y),r,fill=False,edgecolor='#647078',linewidth=.55))

def composition_rows():
    rows=read(OUT/'species_mean_read_composition.csv')
    meta={'species','Status','n_total','n_eligible_ge100'}
    result=[]
    for row in rows:
        vals={g:0. for g in KEYS}
        for k,v in row.items():
            if k not in meta and v:vals[bucket(k)]+=float(v)
        result.append((row,vals))
    return result

def trace_color(color):
    return tuple(.22*c+.78 for c in to_rgb(color))

def detection_pie(ax,x,y,positive,trace,n,color,r=.28):
    if not n:
        ax.add_patch(Circle((x,y),r,facecolor='white',edgecolor='#aaa',hatch='///'))
        ax.text(x+.37,y,'—',va='center',fontsize=8)
        return
    assert positive+trace<=n
    start=90
    for count,fill,hatch in [(positive,color,None),(trace,trace_color(color),'////'),(n-positive-trace,'#e1e4e6',None)]:
        if not count:continue
        end=start+360*count/n
        ax.add_patch(Wedge((x,y),r,start,end,facecolor=fill,edgecolor='#8b9398' if hatch else 'white',linewidth=.35,hatch=hatch))
        start=end
    ax.add_patch(Circle((x,y),r,fill=False,edgecolor='#647078',linewidth=.55))
    ax.text(x+.37,y,f'({positive}+{trace})/{n}',va='center',fontsize=7.7)

def publication_host_tree(tree, rows, tip_rows, yy, xx, depth, composition, stats,
                          traces, selected, observed, mapped_row_count):
    """Render main Figure 2 with essential labels and recorded-country flags."""
    specimens={}
    for r in observed:
        metadata=(r['species'],r['Status'],r['Country'],r['eligible_ge100'])
        assert r['ID'] not in specimens or specimens[r['ID']]==metadata
        specimens[r['ID']]=metadata
    countries=defaultdict(lambda:defaultdict(list))
    for specimen,(species,status,country,eligible) in specimens.items():
        assert country in LOCALITIES, f'Country icon needs a recorded-locality definition: {country}'
        countries[(species,status)][country].append((specimen,eligible=='True'))

    country_audit=[]
    for species,status in rows:
        records=countries[(species,status)]
        assert sum(len(v) for v in records.values())==int(depth[(species,status)]['n_specimens'])
        assert sum(e for v in records.values() for _,e in v)==int(depth[(species,status)]['n_ge100'])
        for country,samples in sorted(records.items()):
            n_eligible=sum(eligible for _,eligible in samples)
            country_audit.append({'species':species,'Status':status,'recorded_country':country,
                'display_code':LOCALITIES[country],'flag_country': 'France' if country=='French Guiana' else
                    ('Federated States of Micronesia' if country=='Micronesia' else country),
                'flag_resolution_source': 'Specimen-origin clarification, 2026-10-08: Micronesia denotes the Federated States of Micronesia.'
                    if country=='Micronesia' else 'Recorded specimen country',
                'n_sampled':len(samples),'n_eligible_ge100':n_eligible,
                'dagger_below_depth_only':n_eligible==0,
                'specimen_IDs':';'.join(s for s,_ in samples),
                'eligible_specimen_IDs':';'.join(s for s,e in samples if e)})

    # Verify that this layout preserves every established plotted composition
    # and detection value, rather than silently recalculating a different rule.
    prior_composition={(r['species'],r['Status']):r for r in read(OUT/'host_tree_composition_pie_data.csv')}
    prior_detection={(r['species'],r['Status'],r['group']):r for r in read(OUT/'host_tree_detection_pie_data.csv')}
    for species,status in rows:
        meta,values=composition[(species,status)]
        old=prior_composition[(species,status)]
        assert int(old['n_eligible'])==int(meta['n_eligible_ge100'])
        assert abs(sum(values.values())-1)<1e-9
        for group,value in values.items():assert abs(float(old[group])-value)<1e-12
        for group in selected:
            key=(species,status,group);r=stats[key];old=prior_detection[key]
            assert int(old['n_eligible'])==int(r['n_eligible_ge100'])
            assert int(old['n_screen_positive'])==int(r['n_positive'])
            assert int(old['n_trace_only'])==traces.get(key,0)
    assert len(prior_composition)==len(rows)
    assert len(prior_detection)==len(rows)*len(selected)

    rowstep=1.12
    fig,ax=plt.subplots(figsize=(22.5,23.1))
    X=lambda value:value*2.1
    Y=lambda value:value*rowstep
    for i,(species,status) in enumerate(rows):
        if species=='A. bruennichi':ax.axhspan(Y(i)-.53,Y(i)+.53,color='#fff5ea',zorder=0)
        elif i%2==0:ax.axhspan(Y(i)-.53,Y(i)+.53,color='#f5f7f8',zorder=0)
    for c in tree.find_clades(order='preorder'):
        if c.clades:
            ax.plot([X(xx[c]),X(xx[c])],[Y(min(yy[k] for k in c.clades)),Y(max(yy[k] for k in c.clades))],color='#667681',lw=1)
            for k in c.clades:ax.plot([X(xx[c]),X(xx[k])],[Y(yy[k]),Y(yy[k])],color='#667681',lw=1)
    for t,indices in tip_rows.items():
        if len(indices)>1:
            ax.plot([X(3.8),X(3.8)],[Y(indices[0]),Y(indices[-1])],color='#667681',lw=1)
            for i in indices:ax.plot([X(3.8),X(3.94)],[Y(i),Y(i)],color='#667681',lw=1)

    group_x=[26.0+3.3*j for j in range(len(selected))]
    for i,(species,status) in enumerate(rows):
        row_y=Y(i);d=depth[(species,status)]
        ax.text(X(4.05),row_y,species,va='center',fontsize=10.7,fontstyle='italic',fontweight='bold' if species=='A. bruennichi' else 'normal')
        ax.text(X(6.2),row_y,status,ha='center',va='center',fontsize=8.7)
        recorded=sorted(countries[(species,status)].items())
        for j,(country,samples) in enumerate(recorded):
            flag_x=17.8+(j-(len(recorded)-1)/2)*1.21
            draw_flag(ax,country,flag_x,row_y,eligible=any(e for _,e in samples))
        ax.text(21.9,row_y,f"{d['n_ge100']}/{d['n_specimens']}",ha='center',va='center',fontsize=9.5)
        meta,values=composition[(species,status)]
        if int(d['n_ge100']):pie(ax,23.65,row_y,list(values.values()),list(COLORS.values()),r=.48)
        else:ax.add_patch(Circle((23.65,row_y),.48,facecolor='white',edgecolor='#aaa',hatch='///'))
        for group,x in zip(selected,group_x):
            key=(species,status,group);r=stats[key]
            detection_pie(ax,x,row_y,int(r['n_positive']),traces.get(key,0),int(r['n_eligible_ge100']),COLORS[group],r=.30)
    if mapped_row_count<len(rows):
        ax.axhline(Y(mapped_row_count-.5),color='#aab8bf',ls='--',lw=1)
        ax.text(X(1.8),Y(mapped_row_count)+.1,'Unplaced species',fontsize=10,ha='center',va='center',color='#586970')
    headings=[(X(1.8),'Host topology'),(X(4.9),'Host species'),(X(6.2),'Sample set'),
        (17.8,'Sampled countries'),(21.9,'Eligible /\nsampled'),(23.65,'Total\nmicrobiome\ncomposition'),
        (group_x[0],'DUSA-like\nOTU2'),(group_x[1],'DUSA relative\nOTU287'),(group_x[2],'Mesoplasma'),
        (group_x[3],'Wolbachia'),(group_x[4],'Spiroplasma')]
    for x,label in headings:ax.text(x,-1.15,label,ha='center',va='center',fontsize=9.8,fontweight='bold')
    ax.set_xlim(-.4,41.45);ax.set_ylim(Y(len(rows)-.4),-2.0);ax.set_aspect('equal');ax.axis('off')
    composition_handles=[Patch(facecolor=color,label=NAMES[group]) for group,color in COLORS.items()]
    fig.legend(handles=composition_handles,loc='lower center',bbox_to_anchor=(.5,.058),ncol=4,frameon=False,fontsize=9.4)
    detection_handles=[Patch(facecolor='#697f91',label='Screen-positive'),
        Patch(facecolor=trace_color('#666666'),edgecolor='#8b9398',hatch='////',label='Trace-only'),
        Patch(facecolor='#e1e4e6',label='No candidate reads'),
        Patch(facecolor='white',edgecolor='#78848b',linestyle='--',label='† Below-depth-only country')]
    fig.legend(handles=detection_handles,loc='lower center',bbox_to_anchor=(.5,.031),ncol=4,frameon=False,fontsize=9.4)
    fig.subplots_adjust(top=.985,bottom=.095,left=.025,right=.995)
    target=ROOT/'figures/main/fig02_host_phylogeny.png'
    target.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(target,dpi=260,facecolor='white');plt.close(fig)
    with (OUT/'host_tree_sample_country_data.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=country_audit[0]);writer.writeheader();writer.writerows(country_audit)
    inputs=['specimen_candidate_screen.csv','species_depth_audit.csv','species_mean_read_composition.csv',
            'species_candidate_prevalence.csv','host_tree_composition_pie_data.csv','host_tree_detection_pie_data.csv',
            'host_phylogeny/host_species_tree_for_plot.nwk','host_phylogeny/species_tip_mapping.csv']
    audit={'figure':'figures/main/fig02_host_phylogeny.png','sample_strata':len(rows),
        'specimens':len(specimens),'mapped_species':len(tip_rows),'unplaced_species':len(rows)-mapped_row_count,
        'composition_cells_verified_unchanged':len(rows)*len(KEYS),
        'detection_pies_verified_unchanged':len(rows)*len(selected),
        'total_microbiome_pie_radius':.48,
        'flags_describe':'All recorded specimen countries; dagger marks countries with no eligible profile.',
        'country_metadata_limits':{'French Guiana':'French national flag; GF locality code.'},
        'country_metadata_resolutions':{'Micronesia':{
            'source_wording':'Micronesia','display_country':'Federated States of Micronesia','display_code':'FM',
            'resolution_source':'Specimen-origin clarification dated 2026-10-08.',
            'resolution_date':'2026-10-08'}},
        'input_sha256':{name:hashlib.sha256((OUT/name).read_bytes()).hexdigest() for name in inputs},
        'output_sha256':hashlib.sha256(target.read_bytes()).hexdigest()}
    (OUT/'main_figure_2_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(target)
    print('MAIN FIGURE 2:',len(rows),'sample strata;',len(specimens),'specimens; all composition and detection values unchanged')


def host_tree():
    treepath=OUT/'host_phylogeny/host_species_tree_for_plot.nwk'
    mappingpath=OUT/'host_phylogeny/species_tip_mapping.csv'
    if not treepath.exists() or not mappingpath.exists():
        raise FileNotFoundError('Missing retained host tree or species mapping')
    tree=Phylo.read(treepath,'newick')
    mapping={r['published_tip']:r['sample_species'] for r in read(mappingpath) if r.get('published_tip') and r.get('mapping_status')=='mapped'}
    depth={(r['species'],r['Status']):r for r in read(OUT/'species_depth_audit.csv')}
    composition={(r['species'],r['Status']):(r,values) for r,values in composition_rows()}
    tips=tree.get_terminals()
    assert all(t.name in mapping for t in tips)
    matched=[mapping[t.name] for t in tips]
    assert len(matched)==len(set(matched))
    stats={(r['species'],r['Status'],r['group']):r for r in read(OUT/'species_candidate_prevalence.csv')}
    observed=read(OUT/'specimen_candidate_screen.csv')
    traces={}
    for r in observed:
        key=(r['species'],r['Status'],r['group'])
        if r['eligible_ge100']=='True' and int(r['curated_group_reads'])>0 and r['screen_positive']=='False':
            traces[key]=traces.get(key,0)+1
    selected=['DUSA-like OTU2','DUSA close sequence relatives','Mesoplasma','Wolbachia','Spiroplasma']
    unknown=sorted({s for s,status in depth if status=='Museum'}-set(matched))
    rows=[];tip_rows={}
    for t in tips:
        species=mapping[t.name]
        indices=[]
        for status in ('Museum','Contemporary'):
            if (species,status) in depth:
                indices.append(len(rows));rows.append((species,status))
        tip_rows[t]=indices
    mapped_row_count=len(rows)
    rows.extend((species,'Museum') for species in unknown)
    yy={t:(idx[0]+idx[-1])/2 for t,idx in tip_rows.items()}
    def y(clade):
        if clade not in yy:yy[clade]=(y(clade.clades[0])+y(clade.clades[-1]))/2
        return yy[clade]
    y(tree.root)
    ranks=tree.depths(unit_branch_lengths=True)
    maxrank=max(ranks[t] for t in tips)
    xx={c:ranks[c]/maxrank*3.8 for c in ranks}
    for t in tips:xx[t]=3.8
    publication_host_tree(tree, rows, tip_rows, yy, xx, depth, composition,
                          stats, traces, selected, observed, mapped_row_count)


if __name__ == '__main__':
    host_tree()
