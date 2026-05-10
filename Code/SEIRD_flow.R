# ============================================================
# SEIRD flow diagram for the spread of HPAI virus in Newfoundland seabirds
# ============================================================

# ============================================================
# Libraries
# ============================================================
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# ============================================================
# SEIRD compartmental structure
# ============================================================
diagram <- grViz("
digraph SEIRD {

  graph [layout = dot, rankdir = LR]

  node [shape = rectangle,
        style = filled,
        fontname = 'Times',
        fontsize = 16,
        fontcolor = black]

  S [label=<S&nbsp;<font point-size='10'><sub>i,j</sub></font>>, fillcolor='#4A90E2']
  E [label=<E&nbsp;<font point-size='10'><sub>i,j</sub></font>>, fillcolor='#F5A623']
  I [label=<I&nbsp;<font point-size='10'><sub>i,j</sub></font>>, fillcolor='#D0021B']
  R [label=<R&nbsp;<font point-size='10'><sub>i,j</sub></font>>, fillcolor='#7ED321']
  D [label=<D&nbsp;<font point-size='10'><sub>i,j</sub></font>>, fillcolor='#9013FE']

  edge [fontname = 'Times',
        fontsize = 14,
        fontcolor = black,
        color = black]

  S -> E [label=<β<font point-size='10'><sub>i,j</sub></font>(t)&nbsp;I&nbsp;<font point-size='10'><sub>i,j</sub></font>&nbsp;/&nbsp;N<font point-size='10'><sub>i,j</sub></font>>]
  E -> I [label=<α>]
  I -> R [label=<γ>]
  I -> D [label=<μ<font point-size='10'><sub>j</sub></font>>]

}
")


diagram

# ============================================================
# Export figure
# ============================================================
svg <- export_svg(diagram)

rsvg_png(
  charToRaw(svg),
  "SEIRD_flow_diagram.png",
  width = 1200,
  height = 800
)
