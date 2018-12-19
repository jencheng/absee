absee
============

[absee](https://rubygems.org/gems/absee) is a Ruby gem that reads ABIF files (DNA sequencing chromatograms)


Documentation
-------------

Documentation is located in the directories.

Example
-------

	%irb
	>> require ‘absee’
	=> true
	>> my_variable = ABSee.new()
	=> #<ABSee:0x000001008599d0>
	>> my_variable.read("/Users/Jenny/Desktop/my_sequence.ab1")
	=> nil
	>> my_variable.get_calledSequence()

Class Methods
-------------

* read(file_location)
	* returns nil
* get_traceA()
	* returns an array with the trace data for adenine
* get_traceG()
	* returns an array with the trace data for guanine
* get_traceC()
	* returns an array with the trace data for cytosine
* get_traceT()
	* returns an array with the trace data for thymine
* get_calledSequence()
	* returns an array with the Basecalled sequence
* get_qualityScores()
	* returns an array with the Basecalled quality scores
* get_peakIndexes()
	* returns an array with indexes of the called sequence in the trace


License
-------

Code is licensed under MIT license 2012.
See the LICENSE file for details.


Credits
-------
Created by Jenny Cheng
