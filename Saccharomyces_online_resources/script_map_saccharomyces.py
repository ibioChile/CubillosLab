import folium
from folium import IFrame
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import base64
import io


df = pd.read_csv("sacc_location.csv")

resolution, width, height = 100, 6, 6

sns.set_style("darkgrid")

def plot_bargraph(df, loc):
    fig, ax = plt.subplots(figsize=(width, height))
    df = df.groupby('Species')['Clade'].value_counts()
    df = df.rename('Samples').reset_index()
    #ax = sns.barplot(x='Species', y='Samples', hue='Clade', hue_order = ["PB1", "PB2","PB3","ADM","S. uvarum"], order = ["S. eubayanus", "S. uvarum"], data=df)
    ax = sns.barplot(x='Clade', y='Samples', order = ["PB1", "PB2","PB3","ADM"], data=df)
    ax.legend().set_visible(False)
    #plt.legend(loc='upper right')
    #ax.set_xlabel("Species")
    ax.set_xlabel("Clade")
    #ax.set_xticklabels(["S. eubayanus", "S. uvarum"], style='italic')
    ax.set_ylabel("Number of Samples")
    ax.set_yticks(np.arange(0, max(df['Samples']) + 1, step=1))
    png = '/Users/pamelacamejo/Documents/IBIO/Francisco Cubillos/Website/plot_{}.png'.format(loc)
    fig.savefig(png, dpi=resolution)

def add_loc_info(loc_info, loc, map):
    loc_info = loc_info.reset_index()
    species = []
    for specie, strain, clade, link in zip(loc_info['Species'],loc_info['Strain'],loc_info['Clade'], loc_info['Link']):
        link_ncbi = '<a href={}>{}</a>'.format(link, strain)
        species.append([specie, clade, link_ncbi])
    table_species = pd.DataFrame(species, columns=["Species","Clade","Strain"])
    table_species['Species'] = [f'<i>{x}</i>' for x in table_species['Species']]

    # Convert table to html
    str_io = io.StringIO()
    table_species.to_html(buf=str_io, classes='table table-striped', border=0.5, col_space=100, justify="center", escape=False, index=False)
    html_str = str_io.getvalue().replace('table border="0.5"', 'table border="0.5" align="center"')
    html_str = html_str.replace('<td>','<td style = "text-align:center">')
    encoded = base64.b64encode(open('/Users/pamelacamejo/Documents/IBIO/Francisco Cubillos/Website/plot_{}.png'.format(loc), 'rb').read()).decode()

    html = f'''<TABLE BORDER=0> \
        <TR>
        <TD>
        <img ALIGN="Right" src="data:image/png;base64,{encoded}" width="400" height="400" > \
        </TD>
        <TD>
        <h1 style="text-align: center; vertical-align: middle;"> Location: {loc} <br/><br/> \
        {html_str}
        </TD>
        </TR>
        </TABLE>'''

    iframe = IFrame(html, width=(width * resolution) + 200, height=(height * resolution) - 100)

    marker = folium.Marker(location=[loc_info['LAT'][0], loc_info['LON'][0]],
                               popup=folium.Popup(iframe, max_width=2650),
                               icon=folium.Icon(icon_color='green',icon="ok"))
    marker.add_to(map)

map = folium.Map(location=[df['LAT'].mean(), df['LON'].mean()], zoom_start=5, tiles='OpenStreetMap')
fg = folium.FeatureGroup(name="Saccharomyces Locations")

for loc in df.LOCATION.unique():
    loc_info = df.loc[df['LOCATION'] == loc]
    plot_bargraph(loc_info, loc)
    add_loc_info(loc_info,loc,map)


map.add_child(folium.LayerControl())

map.save(outfile='map.html')
