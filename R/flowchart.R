########## TOTAL-AMI IV flowchart

options(OutDec = ".")
# Graph
flowchart <- 
  DiagrammeR::grViz("
digraph flowchart {
# graph, node, definitions
  graph [compound = true, nodesep = 1, ranksep = .25, splines=ortho, penwidth = 2]
  node [fontname = Helvetica, shape = rect, fixedsize = false, width = 4, height = 1, style = 'filled',fillcolor = 'beige', penwidth = 2]
  edge [penwidth = 2]

subgraph cluster_Databases {
peripheries=2
graph[style=dashed]
fontname = Helvetica

		SCB [shape = cylinder, fillcolor = 'gray95', label = 'Statistics Sweden \n (SCB)']
    SWEDEHEART[shape = cylinder, fillcolor = 'gray90', label = 'SWEDEHEART \n (RIKS-HIA subregister)']
    SOS [shape = cylinder, fillcolor = 'gray85', label = 'National Board of Health \n and Welfare (SOS) \n (National patient register and \n prescribed drug register)']
label = 'Data sources'
}

rhia [label = '@@1', shape = folder]

subgraph cluster_process {
peripheries=0

	subgraph cluster_Cases {
peripheries=2
graph[style=dashed]
fontname = Helvetica
		node [style=filled];

		scb -> c1 [arrowhead = none]
    c1 -> faulty_cases-> sos 
    sos -> c2 [arrowhead = none]
    c2 -> cases_no_mi 
    cases_no_mi -> c3 [arrowhead = none]
    c3 -> cases_only_first 
    cases_only_first -> c4 [arrowhead = none] 
    c4 -> cases_with_comparators -> study_pop
# empty nodes (connector points)
c1 [height = 0, width = 0, shape=point]
c2 [height = 0, width = 0, shape=point]
c3 [height = 0, width = 0, shape=point]
c4 [height = 0, width = 0, shape=point]


subgraph {
    rank = 'same'
    scb -> scb_invis [style = invis]
}

subgraph {
    rank = 'same'
    c1 -> faulty_cases_rm
}
subgraph {
    rank = 'same'
    c2 -> cases_with_mi
}
subgraph {
    rank = 'same'
    c3 -> cases_not_first
}
subgraph {
    rank = 'same'
    c4 -> cases_no_comparators
}
subgraph {
    rank = 'same'
    empty [style = invis]
    cases_with_comparators ->  empty [style = invis]
}
subgraph {
    rank = 'same'
    empty2 [style = invis]
    study_pop ->  empty2 [style = invis]
}
		label = 'Cases';
}

	subgraph cluster_comparators {
peripheries=2
graph[style=dashed]
fontname = Helvetica
		color=black;
		node [style=filled];
label = 'Comparators'
subgraph{
		scb_comparators -> o1 [arrowhead = none] 
o1 -> faulty_comparators
faulty_comparators -> o2 [arrowhead = none] 
o2-> sos_comparators
sos_comparators -> o3 [arrowhead = none] 
o3 -> o4 [arrowhead = none] 
o4 -> sos_comparators_2
sos_comparators_2 -> study_comparators [arrowhead = none]  


o1 [height = 0, width = 0, shape=point]
o2 [height = 0, width = 0, shape=point]
o3 [height = 0, width = 0, shape=point]
o4 [height = 0, width = 0, shape=point]
}
subgraph {
rank = same
empty3 [style = invis]
scb_comparators -> empty3 [style = invis, width = .1]
}
subgraph {
rank = same
o1 -> faulty_comparators_rm
}
subgraph {
rank = same
 o3 -> comparators_with_ami
}
subgraph {
rank = same
 o4 -> comparators_removed_doubles
}

}


}

 

# subgraph{
# rank = 'same'
#  study_pop -> study_comparators
# }

scb [label = '@@2']
scb_comparators [label = '@@3']
sos [label = '@@4']
faulty_cases [label ='@@5']
faulty_cases_rm [label ='@@6']
faulty_comparators [label = '@@7']
faulty_comparators_rm [label = '@@8']
sos_comparators [label = '@@9']
sos_comparators_2 [label = '@@11']
comparators_with_ami [label = '@@12']
cases_with_mi [label = '@@13']
cases_no_mi [label = '@@14']
cases_only_first [label = '@@15']
cases_with_comparators [label = '@@16']
cases_no_comparators [label = '@@17']
cases_not_first [label = '@@18']
study_pop [label = '@@19']
study_comparators [label = '@@20']
analysis_set [label = '@@21', width = 4, height = 1.5, shape = cylinder, fillcolor = grey90]
scb_invis [style = invis]
invis1 [style = invis, shape = point, height = 0, width = 0]
comparators_removed_doubles [label = '@@22']

rhia -> { scb_comparators scb}

SWEDEHEART -> invis1 [arrowhead = none]
invis1 -> rhia
study_pop -> analysis_set
study_comparators -> analysis_set
SCB -> scb
SCB -> scb_comparators
SOS -> sos
SOS -> sos_comparators

}

[1]: 'All AMI cases identified from \\n SWEDEHEART/RIKS-HIA \\n n = 483,268'
[2]: 'Cases enriched with SCB data  \\n n = 483,268'
[3]: '1:5 Comparator population created from \\n SCB registers \\n n = 1,883,271'
[4]: 'Cases enriched with SOS \\n and RIKS-HIA data \\n n = 490,940 (Doubles introduced)' 
[5]: 'Valid cases \\n  n = 482,366' 
[6]: 'Obviously invalid cases removed \\n  n = 902 (Dead before study start \\n or registry test person)' 
[7]: 'Valid comparators, n = 1,883,253' 
[8]: 'Obviously invalid cases removed  \\n n = 18 (Dead before study start)' 
[9]: 'Comparators enriched with SOS data \\n n = 1,883,253' 
[10]: 'Diagnoses and medication \\n data from SOS' 
[11]: 'Comparators without prior AMI \\n n = 1,821,279'
[12]: 'Comparators found to have prior AMI diagnosis \\n n = 61,973'
[13]: 'Cases with prior MI \\n n = 135,484'
[14]: 'Cases with confirmed MI \\n but no diagnosis of prior MI \\n n = 355,456'
[15]: 'Only first occurrance of \\n every individual, n = 343,921'
[16]: 'Only cases with comparators, n = 335,748'
[17]: 'Cases without comparators, n = 8,173'
[18]: 'Double entries, n = 12,535'
[19]: 'Study population, n = 335,748'
[20]: 'Comparator population matched on \\n study population \\n n = 1,625,734'
[21]: 'Study dataset \\n n = 1,961,144'
[22]: 'Double entries, n = 1'

")
flowchart
flowchart %>% 
  DiagrammeRsvg::export_svg() %>% 
  charToRaw() %>% 
  rsvg::rsvg_pdf("output/flowchart.pdf")
