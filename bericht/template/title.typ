#import "colour.typ": *

#let title_page(
  title: "",
  author: "",
  chair: "",
  school: "",
  degree: "",
  examiner: "",
  supervisor: "",
  submitted: ""
) = {
  page(
    paper: "a4",
    margin: (
      top: 5cm,
      bottom: 3cm,
      x: 2cm,
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
        align(bottom + right, image("resources/TUM_Logo_blau.svg", height: 30%))
      )
    ],
    footer: []
  )[
    #set text(
      font: "TUM Neue Helvetica",
      size: 10pt
    )

    #v(1cm)

    #set align(top + left)
    #text(size: 24pt, [*#title*])

    #v(3cm)

    #text(fill: tum_blue, size: 17pt, [*#author*])

    #v(3cm)

    Bericht zur Ableistung der
    #v(1em)
    *#degree*
    #v(1em)
    an der #school der Technischen Universität München.

    #v(3cm)

    *Prüfer:*\ #examiner
    #v(0em)
    *Betreuer:*\ #supervisor
    #v(0em)
    *Eingereicht am:*\ Munich, #submitted
  ]

  pagebreak()
}
