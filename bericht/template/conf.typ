#import "cover.typ": cover_page
#import "title.typ": title_page
#import "contents.typ": contents_page

#let conf(
  title: "",
  author: "",
  chair: "",
  school: "",
  degree: "",
  examiner: "",
  supervisor: "",
  submitted: "",
  doc
) = {
  cover_page(
    title: title,
    author: author,
    chair: chair,
    school: school
  )

  pagebreak(to: "odd")

  title_page(
    title: title,
    author: author,
    chair: chair,
    school: school,
    degree: degree,
    examiner: examiner,
    supervisor: supervisor,
    submitted: submitted
  )

  set heading(numbering: "1.")
  show heading: it => locate(loc => {
    let levels = counter(heading).at(loc)

    set text(font: "TUM Neue Helvetica")
    if it.level == 1 [
      #if levels.at(0) != 1 {
        pagebreak(to: "odd")
      }
      #set text(size: 24pt)
      #if levels.at(0) != 0 {
        numbering("1", levels.at(0))
      }
      #it.body
      #v(1em, weak: true)
    ] else if it.level == 2 [
      #set text(size: 16pt)
      #v(1em)
      #numbering("1.1", levels.at(0), levels.at(1))
      #it.body
      #v(1em, weak: true)
    ] else if it.level == 3 [
      #set text(size: 16pt)
      #v(1em, weak: true)
      #numbering("1.1.1", levels.at(0), levels.at(1), levels.at(2))
      #it.body
      #v(1em, weak: true)
    ] else [
      #set text(size: 12pt)
      #v(1em, weak: true)
      #it.body
      #v(1em, weak: true)
    ]
  })

  set page(
    paper: "a4",
    margin: (
      top: 3cm,
      bottom: 3cm,
      x: 2cm,
    ),
    header: [],
    footer: []
  )

  contents_page()

  pagebreak(to: "odd")

  set par(justify: true)
  set align(left)
  set text(
    font: "Times New Roman",
    size: 12pt
  )

  set page(
    header: [],
    footer: locate(loc =>
      if calc.rem(loc.page(), 2) == 0 {
        align(left, text(font: "TUM Neue Helvetica", size: 10pt, counter(page).display("1")));
      } else {
        align(right, text(font: "TUM Neue Helvetica", size: 10pt, counter(page).display("1")));
      }
    )
  )

  doc
}

