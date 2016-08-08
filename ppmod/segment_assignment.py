"""Module that contains a Jupyter interactive GUI for assigning segments 
to a topology"""
from __future__ import print_function, division, absolute_import, unicode_literals

import ipywidgets as widgets
import ppmod.topology as t
from collections import OrderedDict, namedtuple

parallel_pairs_list_def = ["P1:P2", "P3:P4", "P5:P6","P7:P8","P9:P10","P11:P12", "GCNsh:GCNsh"]
antiparallel_pairs_list_def =["APHsh:APHsh","BCR:BCR", "APH4:APH4"]
#TODO: APH4 is just a guess
segment_strengths_parallel_def = ["P3:P4", "P1:P2", "P9:P10", "P11:P12", "P5:P6", "GCNsh:GCNsh", "P7:P8"]
segment_strengths_antiparallel_def = ["APHsh:APHsh", "BCR:BCR", "APH4:APH4"]
types_list_def = ['SN', 'S', 'A']



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


def gui(topology, parallel_pairs_list=parallel_pairs_list_def,
                  antiparallel_pairs_list=antiparallel_pairs_list_def,
                  segment_strengths_parallel=segment_strengths_parallel_def,
                  segment_strengths_antiparallel=segment_strengths_antiparallel_def,
                  types_list=types_list_def,
                  auto_display=True,
                  auto_assign=True):
    """Returns a displayable Jupyer notebook GUI for for assigning segments 
       to a topology and a has dictionary with the most important widgets"""                
    
    #initialization
    pairs = t.get_complete_pairs_dict_from_topology(topology)
    
    
    
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
            ass_dict = t.segment_assignments_to_dict(asignment_textbox.value)        
        except ValueError:
            print("Warning: some segments not assigned!")   
            
            
    ###########################################################################
    #MAIN BODY
    ###########################################################################
    
    #Create dropdowns
    for p, pp in pairs.iteritems():
        if is_pair_parallel(pp):
            option_list = ['None'] + parallel_pairs_list
        else:        
            option_list = ['None'] + antiparallel_pairs_list
    
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
        widgets.VBox(display_widget_list),
        automatic_button,
        save_results_button,
        asignment_textbox,
    
            
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