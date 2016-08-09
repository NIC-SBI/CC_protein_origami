"""Module that contains a Jupyter interactive GUI for assigning segments 
to a topology"""
from __future__ import print_function, division, absolute_import, unicode_literals

import ipywidgets as widgets
import ppmod.topology as t
from collections import OrderedDict, namedtuple
import ppmod.utils as u
import collections

pairs_parallel_def = ["P1:P2", "P3:P4", "P5:P6","P7:P8","P9:P10","P11:P12", "GCNsh:GCNsh"]
pairs_antiparallel_def =["APHsh:APHsh","BCR:BCR", "APH4:APH4"]
#TODO: APH4 is just a guess
segment_strengths_parallel_def = ["P3:P4", "P1:P2", "P9:P10", "P11:P12", "P5:P6", "GCNsh:GCNsh", "P7:P8"]
segment_strengths_antiparallel_def = ["APHsh:APHsh", "BCR:BCR", "APH4:APH4"]
types_list_def = ['SN', 'S', 'A']

def get_pairs_from_topology(topology):
    """Returns the pairs A, B, C... from the topology string or list"""
    tU = [p.upper() for p in topology]    
        
    return sorted(list(set(list(tU))))
    
def get_complete_pairs_dict_from_topology(topology):
    """Returns the dictionary A-> (A,a), B->(B, B) ... from the topology string or list"""
    pair_keys = collections.OrderedDict()    
        
    for p in topology:
        P=p.upper()
        l = pair_keys.get(P, [])
        l.append(p)
        pair_keys[P] = l
        
    return pair_keys  
    
    

def segment_assignments_to_dict(rep_str):
    """Parses segment_assignments rules into a dictionary. The string is in the form

    A->SEG1:SEG2
    B->SEG:SEG    
    """
    rep_str = rep_str.strip(" \n")
    rep_lines = rep_str.split("\n")
    reps = collections.OrderedDict()
    for line in rep_lines:     
        if not "->" in line: 
            #print("skipping line", line)            
            continue        
        k,v = line.split("->")
        k = k.strip().upper()
        v = v.strip()
       
        v1,v2 = v.split(":")
        v1 = v1.strip()
        v2 = v2.strip()
        reps[k] = [v1,v2]
    return reps

def do_assignment_replacements(topology, assignments):
    """Replaces topology segments with replacemnt rules. 
       The first occurance of a name is replaced witht he first element in the list"""
    if u.is_str(assignments):
        assignments = segment_assignments_to_dict(assignments) 
    else:
        assignments = assignments.copy()
    
    keys = assignments.keys()    
    ret =[]
    for s in topology:
        sl = s.upper()
        if sl in keys:    
            ret.append(assignments[sl].pop(0))
        else:
            ret.append(s)
    return ret    

def seq_to_seq_map(dict_or_file_or_dataframe):
    """Takes a dictionary (in which case it does not modify it) or a dataframe, 
    or an excel file name and returns a dict mapping segment names to sequences.
    Columns must be PID and Sequence"""
    import pandas as pd
    #if is dictionary just return it
    if isinstance(dict_or_file_or_dataframe, dict):
        return dict_or_file_or_dataframe
        
    #if is file, load it to a dataframe
    if u.is_str(dict_or_file_or_dataframe):     
        df = pd.read_excel(dict_or_file_or_dataframe)
    else:
        df = dict_or_file_or_dataframe
    #else it must be a dataframe. Do your magic:)
    df=df.dropna()
    
    #cleanup whitespace
    df.PID=df.PID.str.strip()
    
    peps = {}
    for n, pep in df.iterrows():

        peps[pep.PID] = pep.Sequence.replace('-', '').strip()
    return peps
    
    
def get_annotated_sequence(segments, seg_to_seq, linkers="SGPGS", N_tag="", C_tag=""):
    """Returns an anottated amino acid sequnce of the poylhedra.
    Linkers can be a single string or an array of string with the required number of linkers.
    N_tag, C_tag - are appended to the left and right side of the string"""

    seg_to_seq = seq_to_seq_map(seg_to_seq)
    N = len(segments)    
    if u.is_str(linkers):
        linkers = [linkers]*(N-1)
        
    assert len(linkers)==N-1, ("Length of linkers must be one less than the number of segments."+
                               "Is {NL}, but should be {N}".format(NL=len(linkers), N=N))
                               
    aa_segments = [seg_to_seq[s]  +"\t|"+s for s in segments]
    
    lines = [N_tag] + list(u.roundrobin(aa_segments, linkers)) + [C_tag]
    lines = "\n".join(lines)
    return lines
    
def deannotate_sequence(annotated_seq, remove_whitespace=False):
    """Deletes the annotation of segments. Optionally removes whitespace.
    Either a tab or a | are sufficent, both are also acceptable"""
    import re    
    #mathces from a tab or | or both to the end of line.    
    reg_exp = re.compile(r"\s*(\t|\||\t\s*\|).*\n", re.MULTILINE) 
    res = re.sub(reg_exp,"\n",annotated_seq)    
    
    if remove_whitespace:
        res=res.replace(" ","")
        res=res.replace("\n","")
        res=res.replace("\t","")
    return res
    
def is_pair_parallel(pair):
    return pair[0]==pair[1]

def pair_description(pair):
    return str(pair[0])+':'+str(pair[1])

    
def splice_in_type(pair, type):
    """Takes a pair and splices in the type APH:APH -> APHshSN:APHshSN"""
    sp = pair.split(':')
    if len(sp) == 2:
        p1,p2=pair.split(':')
        return p1+type+':'+p2+type
    else:
        return pair

GUI_result = namedtuple("GUI_result", "gui pairs_dict pair_dropdowns, type_dropdowns, result_text") 

#TODO
# - Document
# - write TCO at the end
# - add linker input box

def segment_assignment_gui(topology, pairs_parallel=pairs_parallel_def,
                  pairs_antiparallel=pairs_antiparallel_def,
                  segment_strengths_parallel=segment_strengths_parallel_def,
                  segment_strengths_antiparallel=segment_strengths_antiparallel_def,
                  types_list=types_list_def,
                  auto_display=True,
                  auto_assign=True):
    """Returns a displayable Jupyer notebook GUI for for assigning segments 
       to a topology and a has dictionary with the most important widgets"""                
    
    #initialization
    pairs = get_complete_pairs_dict_from_topology(topology)
    
    
    
    widget_dropdown_dict = OrderedDict()
    widget_type_dict = OrderedDict()
    display_widget_list = []    
    
    #helper events and functions  
    def on_value_change(change):
        print(change['owner'].view_name)    
    
    def automatic_assignment_click(btn):
        distances = t.get_segment_distances_dict(topology)
        distances = OrderedDict(sorted(distances.items(), key=lambda t: t[1], reverse=True))
        
        
        p_counter, a_counter = 0, 0
        
        for p in distances.keys():
            if is_pair_parallel(pairs[p]):
                widget_dropdown_dict[p].value=segment_strengths_parallel[p_counter]
                p_counter += 1
            else:
                widget_dropdown_dict[p].value=segment_strengths_antiparallel[a_counter]
                a_counter += 1
    
        save_results_click(btn)        
    
    def save_results_click(btn):   
        lines  = []
        for p, pp in pairs.iteritems():
            pair = widget_dropdown_dict[p].value
            pair = splice_in_type(pair, widget_type_dict[p].value)
            l = "{p}->{pair}".format(p=p, pair=pair)
            lines.append(l)
        asignment_textbox.value="\n".join(lines)
    
        
        if widget_dropdown_dict[topology[0].upper()].value in segment_strengths_parallel[-2:]:
            print("Warning: first segment has weak binding! It is not advisable to put it at the beginning of the chain.")
        if widget_dropdown_dict[topology[-1].upper()].value in segment_strengths_parallel[-2:]:
            print("Warning: las segment has weak binding! It is not advisable to put it at the end of the chain.")              
                  
        try:
            segment_assignments_to_dict(asignment_textbox.value)        
        except ValueError:
            print("Warning: some segments not assigned!")   
            
            
    ###########################################################################
    #MAIN BODY
    ###########################################################################
    
    #Create dropdowns
    for p, pp in pairs.iteritems():
        if is_pair_parallel(pp):
            option_list = ['None'] + pairs_parallel
        else:        
            option_list = ['None'] + pairs_antiparallel
    
        pair_dropdown = widgets.Dropdown(
                        options=option_list,
                        value='None',
                        description=pair_description(pp))
        #pair_dropdown.view_name=p
        widget_dropdown_dict[p]=pair_dropdown
        
        type_dropdown = widgets.Dropdown(options=types_list)
        widget_type_dict[p]=type_dropdown                
                        
        
        display_widget_list.append(widgets.HBox([pair_dropdown, type_dropdown]))
    
    save_results_button = widgets.Button(description='Print assignment')
    save_results_button.on_click(save_results_click)
    
    asignment_textbox = widgets.Textarea()
    
    automatic_button = widgets.Button(description="Automatic asignment")   
    automatic_button.on_click(automatic_assignment_click)
    
    
    gui_list = widgets.VBox([
        widgets.HTML("Topology: "+topology),        
        widgets.VBox(display_widget_list),
        widgets.HBox([automatic_button, save_results_button, ]),
        asignment_textbox    
        ])
    
    result = GUI_result(gui=gui_list, pairs_dict = pairs, 
                        pair_dropdowns=widget_dropdown_dict, 
                        type_dropdowns=widget_type_dict, 
                        result_text=asignment_textbox)
    
    if auto_assign:
        automatic_assignment_click(automatic_button)

    if auto_display:
        from IPython.display import display
        display(gui_list)
        
        
    return result
   
#GUI_text_result = namedtuple("GUI_text_result", "gui result_text") 
   
def text_edit_gui(text="",caption="", auto_display=True):
    """IPython widget GUI to edit a block of text. Returns the text area wiget""" 
        
        
    textbox = widgets.Textarea(text)
    caption = widgets.HTML(caption)    
    gui = widgets.VBox([caption, textbox])    
    #V 5 of ipywigets is needed for this    
    #gui.layout.width  = '100%'    
    #gui.layout.height = '200px'
    if auto_display:
        from IPython.display import display
        display(gui)
    return textbox
    
#TODO test
def get_included_pairs_info(excel_name, sheetname='pairs', included_pairs=None):
    import pandas as pd

    all_pairs = pd.read_excel(excel_name, sheetname=sheetname).dropna(how='all')
        
    
    pairs_info = all_pairs[all_pairs.pair.isin(included_pairs)]
    pairs_info = pairs_info.sort_values(by='strength', ascending=False)
    

    pairs_dict = [OrderedDict(row) for i, row in pairs_info.iterrows()]
    import yaml
    
    #taken from http://stackoverflow.com/a/16782282/952600
    #Support  for pretty printnig of OrderedDict    
    def represent_ordereddict(dumper, data):
        value = []
    
        for item_key, item_value in data.items():
            node_key = dumper.represent_data(item_key)
            node_value = dumper.represent_data(item_value)
    
            value.append((node_key, node_value))
    
        return yaml.nodes.MappingNode(u'tag:yaml.org,2002:map', value)
    
    yaml.add_representer(OrderedDict, represent_ordereddict)
    return yaml.dump(pairs_dict, encoding=None).replace("!!python/unicode ","")
    
make_config_template ='''\
model_name = "{model_name}"
annotated_sequence = """
{annotated_sequence}
"""

pairs_info = """
{pairs_info}     
"""   

if __name__ == "__main__":
    import ppmod.make_json as mj   
    import ppmod.segment_assignment as sa
    entire_sequence = sa.deannotate_sequence(annotated_sequence, remove_whitespace=True)  
    pairs = mj.load_pairs(pairs_info)    
    mj.generate_json(model_name, entire_sequence, annotated_sequence, pairs)

'''

#TODO test
def write_make_config(model_name, annotated_sequence, pairs_info, 
                      out_name='make_config.py'):
    """Writes the make_config file. If outname is None return the string"""
    data = make_config_template.format(model_name=model_name, 
    annotated_sequence=annotated_sequence,
    pairs_info=pairs_info)
    
    if out_name is None:
        return data
    else:
        with open(out_name, "w") as text_file:
            text_file.write(data)