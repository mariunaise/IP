#import "@preview/cetz:0.2.2"
#import "@preview/fletcher:0.5.1"
#import "@preview/gentle-clues:0.9.0"
#import "@preview/glossarium:0.4.1": *
#import "@preview/lovelace:0.3.0"
#import "@preview/tablex:0.0.8"
#import "@preview/unify:0.6.0"
#import "@preview/quill:0.3.0"
#import "@preview/equate:0.2.0": equate
#import "@preview/drafting:0.2.0": *


#show: equate.with(breakable: true, sub-numbering: true)
#set math.equation(numbering: "(1.1)")

#show figure.where(
  kind: table
): set figure.caption(position: top)

#import "template/conf.typ": conf

#show: make-glossary

#set document(title: "Sehr seriöser IP Titel", author: "Marius Drechsler")

#set text(lang: "de")

#show: doc => conf(
  title: "Sehr seriöser IP Titel",
  author: "Marius Drechsler",
  chair: "Lehrstuhl für Sicherheit in der Informationstechnik",
  school: "School of Computation, Information and Technology",
  degree: "Ingenieurspraxis",
  examiner: "Dr. Michael Pehl",
  supervisor: "M.Sc. Jonas Ruchti",
  submitted: "2025",
  doc
)
#include "content/introduction.typ"
#include "content/background.typ"
#include "content/execution.typ"
#include "content/results.typ"

#include "glossary.typ"

#counter(heading).update(0)
#bibliography("bibliography.bib", style: "ieee")

