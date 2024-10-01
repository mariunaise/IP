#import "colour.typ": *

#let cover_page(
  title: "",
  author: "",
  chair: "",
  school: ""
) = {
  page(
    paper: "a4",
    margin: (
      top: 3cm,
      bottom: 1cm,
      x: 1cm,
    ),
    header: [
      #grid(
        columns: (1fr, 1fr),
        rows: (auto),
        text(
          fill: tum_blue,
          size: 8pt,
          font: "TUM Neue Helvetica",
          [#chair \ #school \ Technische Universität München]
        ),
        align(bottom + right, image("resources/TUM_Logo_blau.svg", height: 50%))
      )
    ],
    footer: []
  )[
    #v(1cm)

    #align(top + left)[#text(font: "TUM Neue Helvetica", size: 24pt, [*#title*])]
    
    #v(3cm)

    #text(font: "TUM Neue Helvetica", fill: tum_blue, size: 17pt, [*#author*])
    
    #align(bottom + right)[#image("resources/TUM_Tower.png", width: 60%)]
  ]

  pagebreak()
}
