Overview 
============

Genbank Browser is a tool for the graphic exploration of genbank files. 

It can automatically read a generic genbank file and expose its content as a web-based genome browser, in which tracks and annotation are automatically handled. Features and their annotations are automatically indexed.

This software is written on top of BioPython, Sphinx, and Genome Viewer tools. 


External Dependencies
======================
Most dependencies are included in this repository to ensure compatibility:

 -  Sphinx **(Require compiling)**: cd sphinx/ &&  ./configure && make
 -  Cherrypy 
 -  Bottle 
 -  BCBio  

In addition, you will need Python 2.6+, BioPython and modern web browser (Chrome recommended).


Features
============

- Automatic indexing of fearures and descriptions
- Instant location of genes and features
- Transparent integration with Genome View
- Web service ready: start gbank-browser as a daemon and serve genbank files as a genome browser
- Combine genbank annotations with your custom cross-linking gene and protein ids. 

Roadmap
==============================
- Combine genbank with custom GFF3 annotations
- Easy edition of genbank features 

Test example
===============
  $ python genbank_browser.py --genbank examples/hs_ref_GRCh38_chr3.gbs  --gname gene --refresh start


Online example
==================

http://ct.bork.embl.de/ctbrowser