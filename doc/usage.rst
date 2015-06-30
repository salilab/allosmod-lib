Usage
*****

The AllosMod protocol can be run by means of a command line tool
(``allosmod``). Each component of the protocol is also a Python package,
which can be called directly from other Python software
(via ``import allosmod``).

Overview
========

Running the basic protocol consists of three steps:

#. Create a set of input files specifying the sequence to model, any
   structures to use as templates, and AllosMod parameters.

#. Run ``allosmod setup`` to verify the inputs and generate a script file.

#. Run the script file (typically on a Linux cluster) to set up the AllosMod
   energy landscape and generate Modeller input files to sample it. In some
   cases this sampling is then carried out; in others you will need to run
   Modeller on the input files.
