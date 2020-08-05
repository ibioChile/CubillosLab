import pandas as pd
import numpy as np

from bokeh.layouts import column, layout
from bokeh.models import ColumnDataSource, CustomJS, Select, Div
from bokeh.plotting import Figure, save

# Import data
df = pd.read_csv("/Users/pamelacamejo/Documents/IBIO/Francisco_Cubillos/Website/Bokeh/metadata.csv")
df["color"] = np.where(df["country"] == "Chile", "orange", "blue")

desc = Div(text=open("/Users/pamelacamejo/Documents/IBIO/Francisco_Cubillos/Website/Bokeh/description.html").read(),
           sizing_mode="stretch_width")

# Create Column Data Source that will be used by the plot
original_source = ColumnDataSource(data=df)
source = ColumnDataSource(
    data=dict(strain=[], accession=[], color=[], CO2_lost_gL=[], country=[], location=[], clade=[]))

### Main plot

# Create plot
TOOLTIPS = [
    ("Strain", "@strain"),
    ("Accession", "@accession"),
    ("Country", "@country"),
    ("Location", "@location"),
    ("Clade", "@clade"),
    ("CO\u2082 lost (g/L)", "@CO2_lost_gL")
]

p = Figure(x_range=list(df.strain), plot_height=600, plot_width=1000,
           title="Fermentation rate", tooltips=TOOLTIPS, toolbar_location=None)
p.circle(x="strain", y="CO2_lost_gL", source=source, size=7, color="color", line_color=None, fill_alpha=0.5)
p.xaxis.axis_label = "Strain"
p.yaxis.axis_label = "CO\u2082 lost (g/L)"
p.axis.axis_label_text_font_style = "bold"
p.grid.grid_line_alpha = 0.3
p.xaxis.major_label_orientation = "vertical"

# Select country
available_countries = list(set(df['country']))
available_countries.sort()
# Select location
available_loc = list(set(df['location']))
available_loc.sort()
# Select clade
clades_list = list(set(df['clade']))
available_clade = [x for x in clades_list if str(x) != 'nan']
available_clade.sort()
# Select strain
available_strain = list(set(df['strain']))
available_strain.sort()

code = """
var data = original_source.data;
var source_data = source.data;

var country_data = data['country'];
var location_data = data['location'];
var accession_data = data['accession'];
var strain_data = data['strain'];
var clade_data = data['clade'];
var color_data = data['color'];
var CO2_lost_gL_data = data['CO2_lost_gL']

var selected_countries = cb_obj.value;

var source_country = source_data['country'];
source_country.length = 0;

var source_location = source_data['location'];
source_location.length = 0;

var source_strain = source_data['strain'];
source_strain.length = 0;

var source_clade = source_data['clade'];
source_clade.length = 0;

var source_color = source_data['color'];
source_color.length = 0;

var source_CO2_lost_gL = source_data['CO2_lost_gL'];
source_CO2_lost_gL.length = 0;

var source_accession = source_data['accession'];
source_accession.length = 0;

var unique_sets = [];

if(selected_countries != "All"){
        for(var i = 0; i < strain_data.length; i++) {
            if(selected_countries.indexOf(country_data[i]) >= 0){
                source_country.push(country_data[i]);
                source_location.push(location_data[i]);
                source_strain.push(strain_data[i]);
                source_clade.push(clade_data[i]);
                source_CO2_lost_gL.push(CO2_lost_gL_data[i]);
                source_color.push(color_data[i]);
                source_accession.push(accession_data[i]);
                
                if ( !unique_sets.includes(location_data[i]) ) {
                unique_sets.push(location_data[i]);
                }
                }
        }
        source.change.emit();
        unique_sets.push("All")
        location.options = unique_sets;
        }
else { 
    for(var i = 0; i < strain_data.length; i++) {
            source_country.push(country_data[i]);
            source_location.push(location_data[i]);
            source_strain.push(strain_data[i]);
            source_clade.push(clade_data[i]);
            source_CO2_lost_gL.push(CO2_lost_gL_data[i]);
            source_color.push(color_data[i]);
            source_accession.push(accession_data[i]);
            }

    source.change.emit();
    var options_all = [...new Set(location_data)];
    options_all.push("All")
    location.options = options_all
    }
"""
code2 = """
var data = original_source.data;
var source_data = source.data;

var country_data = data['country'];
var location_data = data['location'];
var strain_data = data['strain'];
var clade_data = data['clade'];
var color_data = data['color'];
var CO2_lost_gL_data = data['CO2_lost_gL'];
var accession_data = data['accession'];

var selected_location = cb_obj.value;

var source_country = source_data['country'];
source_country.length = 0;

var source_location = source_data['location'];
source_location.length = 0;

var source_strain = source_data['strain'];
source_strain.length = 0;

var source_clade = source_data['clade'];
source_clade.length = 0;

var source_color = source_data['color'];
source_color.length = 0;

var source_CO2_lost_gL = source_data['CO2_lost_gL'];
source_CO2_lost_gL.length = 0;

var source_accession = source_data['accession'];
source_accession.length = 0;

var unique_sets= []

if(selected_location != "All"){
        for(var i = 0; i < strain_data.length; i++) {
            if(selected_location.indexOf(location_data[i]) >= 0){
                source_country.push(country_data[i]);
                source_location.push(location_data[i]);
                source_strain.push(strain_data[i]);
                source_clade.push(clade_data[i]);
                source_CO2_lost_gL.push(CO2_lost_gL_data[i]);
                source_color.push(color_data[i]);
                source_accession.push(accession_data[i]);

                if ( !unique_sets.includes(clade_data[i]) ) {
                unique_sets.push(clade_data[i]);
                }
            }
        }
        source.change.emit();
        unique_sets.push("All")
        clade.options = unique_sets;
        }
else { 
        for(var i = 0; i < strain_data.length; i++) {
            source_country.push(country_data[i]);
            source_location.push(location_data[i]);
            source_strain.push(strain_data[i]);
            source_clade.push(clade_data[i]);
            source_CO2_lost_gL.push(CO2_lost_gL_data[i]);
            source_color.push(color_data[i]);
            source_accession.push(accession_data[i]);
            }
        source.change.emit();
        var options_all = [...new Set(clade_data)];
        options_all.push("All")
        clade.options = options_all
    }
"""
code3 = """
var data = original_source.data;
var source_data = source.data;
var loc_val = location.value
console.log(loc_val)
var country_data = data['country'];
var location_data = data['location'];
var strain_data = data['strain'];
var clade_data = data['clade'];
var color_data = data['color'];
var CO2_lost_gL_data = data['CO2_lost_gL'];
var accession_data = data['accession'];

var selected_clade = cb_obj.value;

var source_country = source_data['country'];
source_country.length = 0;

var source_location = source_data['location'];
source_location.length = 0;

var source_strain = source_data['strain'];
source_strain.length = 0;

var source_clade = source_data['clade'];
source_clade.length = 0;

var source_color = source_data['color'];
source_color.length = 0;

var source_CO2_lost_gL = source_data['CO2_lost_gL'];
source_CO2_lost_gL.length = 0;

var source_accession = source_data['accession'];
source_accession.length = 0;

var unique_sets= []

if(selected_clade != "All"){
        for(var i = 0; i < strain_data.length; i++) {
            if(selected_clade.indexOf(clade_data[i]) >= 0 & loc_val.indexOf(location_data[i]) >= 0 ){
                source_country.push(country_data[i]);
                source_location.push(location_data[i]);
                source_strain.push(strain_data[i]);
                source_clade.push(clade_data[i]);
                source_CO2_lost_gL.push(CO2_lost_gL_data[i]);
                source_color.push(color_data[i]);
                source_accession.push(accession_data[i]);
                
                if ( !unique_sets.includes(strain_data[i]) ) {
                unique_sets.push(strain_data[i]);
                }
            }
        }
        source.change.emit();
        unique_sets.push("All")
        strain.options = unique_sets;
        }
else { 
        for(var i = 0; i < strain_data.length; i++) {
            source_country.push(country_data[i]);
            source_location.push(location_data[i]);
            source_strain.push(strain_data[i]);
            source_clade.push(clade_data[i]);
            source_CO2_lost_gL.push(CO2_lost_gL_data[i]);
            source_color.push(color_data[i]);
            source_accession.push(accession_data[i]);
            }
        source.change.emit();
        var options_all = [...new Set(strain_data)];
        options_all.push("All")
        strain.options = options_all
    }
"""
code4 = """
var data = original_source.data;
var source_data = source.data;

var country_data = data['country'];
var location_data = data['location'];
var strain_data = data['strain'];
var clade_data = data['clade'];
var color_data = data['color'];
var CO2_lost_gL_data = data['CO2_lost_gL'];
var accession_data = data['accession'];

var selected_strain = cb_obj.value;

var source_country = source_data['country'];
source_country.length = 0;

var source_location = source_data['location'];
source_location.length = 0;

var source_strain = source_data['strain'];
source_strain.length = 0;

var source_clade = source_data['clade'];
source_clade.length = 0;

var source_color = source_data['color'];
source_color.length = 0;

var source_CO2_lost_gL = source_data['CO2_lost_gL'];
source_CO2_lost_gL.length = 0;

var source_accession = source_data['accession'];
source_accession.length = 0;

if(selected_strain != "All"){
        for(var i = 0; i < strain_data.length; i++) {
            if(selected_strain.indexOf(strain_data[i]) >= 0){
                source_country.push(country_data[i]);
                source_location.push(location_data[i]);
                source_strain.push(strain_data[i]);
                source_clade.push(clade_data[i]);
                source_CO2_lost_gL.push(CO2_lost_gL_data[i]);
                source_color.push(color_data[i]);
                source_accession.push(accession_data[i]);
                }
        }
        source.change.emit()
        }
else { 
        for(var i = 0; i < strain_data.length; i++) {
            source_country.push(country_data[i]);
            source_location.push(location_data[i]);
            source_strain.push(strain_data[i]);
            source_clade.push(clade_data[i]);
            source_CO2_lost_gL.push(CO2_lost_gL_data[i]);
            source_color.push(color_data[i]);
            source_accession.push(accession_data[i]);
            }
        source.change.emit();
    }
"""

country = Select(title="Country", value="", options=available_countries + ["All"])
location = Select(title="Location", value="", options=available_loc + ["All"])
clade = Select(title="Clade", value="", options=available_clade + ["All"])
strain = Select(title="Strain", value="", options=available_strain + ["All"])

country_callback = CustomJS(args={'source': source, 'original_source': original_source, 'location': location}, code=code)
location_callback = CustomJS(args={'source': source, 'original_source': original_source, 'clade': clade}, code=code2)
clade_callback = CustomJS(args={'source': source, 'original_source': original_source, 'strain': strain,
                                'location': location}, code=code3)
strain_callback = CustomJS(args={'source': source, 'original_source': original_source}, code=code4)

strain.js_on_change('value', strain_callback)
clade.js_on_change('value', clade_callback)
location.js_on_change('value', location_callback)
country.js_on_change('value', country_callback)

controls = [country, location, clade, strain]
inputs = column(*controls, width=320, height=1000)
inputs.sizing_mode = "fixed"
l = layout([
    [desc],
    [inputs, p],
])

# Export plot
save(l, filename="/Users/pamelacamejo/Documents/IBIO/Francisco_Cubillos/Website/Bokeh/Phenotype_sebuyanus.html",
     title="S.eubayanus Phenotype")
