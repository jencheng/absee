# absee
#
# Jenny Cheng
# jencheng@ginkgobioworks.com
#
# based off of Abi.cs by Ronaldo Rodrigues Ferreira
#
# extracts the data from ABIF files
#
# MIT license 2012

class ABSee

  #variables
  @traceA = []
  @traceG = []
  @traceC = []
  @traceT = []
  @calledSequences = []
  @peakIndexes = []
  @qualityScores = []
  
  #opens the ABIF sequencing / chromatogram file
  #checks for ABIF file type
  #major ABIF versions greater than 1 are not supported
  #
  #== Parameters:
  #filename::
  #   a string containing the filename (including the path and extensions)
  #
  #== Returns:
  #  Six arrays: trace data for A, C, G, T, called sequence, and peak indexes
  def read(filename)
    if filename.instance_of?(StringIO) or filename.instance_of?(File)
      abFile = filename
    elsif filename.instance_of? String
      #opens ab1 as a File object
      abFile = open(filename)
    else
      raise "Don't know what to do with input type #{filename.class}"
    end

    byteArray = ""
    #// here we read the first four bytes. It is important
    #// to remember that we do not seek back the file, just
    #// because it is not necessary to do this.
    abFile.seek(0, IO::SEEK_SET)
    abFile.read(4, byteArray)
    #ABIF file indicator
    if byteArray == "ABIF"
      processAB(abFile)
    else
      raise "file not recognized as ABIF"
    end
  end


  ##accessors
  #== Returns:
  #  an array with the trace data for adenine
  def get_traceA()
    return @traceA
  end
  #== Returns:
  #  an array with the trace data for guanine
  def get_traceG()
    return @traceG
  end
  #== Returns:
  #  an array with the trace data for thymine
  def get_traceT()
    return @traceT
  end
  #== Returns:
  #  an array with the trace data for cytosine
  def get_traceC()
    return @traceC
  end
  #== Returns:
  #  an array with the Basecalled sequence
  def get_calledSequence()
    return @calledSequence
  end
  #== Returns:
  #  an array with the Basecalled quality scores
  def get_qualityScores()
    return @qualityScores
  end
  #== Returns:
  #  an array with the peak indexes
  def get_peakIndexes()
    return @peakIndexes
  end
  
  
  private
  
  #process the opened ABIF filestream, and calls subsequent methods to extract the data
  #
  #== Parameters:
  #filestream:: an opened File
  #
  #== Returns:
  #Six arrays: trace data for A, C, G, T, called sequence, and peak indexes
  #readAB returns the results of this method
  def processAB(filestream)
    #// here, we can read the ABIF header information
    version = readUnsignedByte_2(4, filestream)
    #// major versions greater than 1 are not supported
    #// Applied Biosystems rules
    if (version / 100 > 1)
      raise "ABIF version #{version} not supported (only supported for version less than 1)"
    end
    #// we just read ABIF, so we don't need more information than that
    numElements = readUnsignedByte_4(18, filestream)
    dataOffset = readUnsignedByte_4(26, filestream)
    directory = readDirectoryEntry(filestream, dataOffset, numElements)
    numSamples, numBases = gatherInformation(directory, numElements)
    samples_a, samples_c, samples_g, samples_t = getSamples(filestream, directory, numElements, numSamples)
    called_sequence = getCalledSequence(filestream, directory, numElements, numBases)
    quality_scores = getQualityScores(filestream, directory, numElements, numBases)
    peak_indexes = getPeakIndexes(filestream, directory, numElements, numBases)
    ##return samples_a, samples_c, samples_g, samples_t, called_sequence, peak_indexes, quality_scores
    @traceA = samples_a
    @traceC = samples_c
    @traceG = samples_g
    @traceT = samples_t
    @calledSequence = called_sequence
    @qualityScores = quality_scores
    @peakIndexes = peak_indexes
    nil
  end
  
  #reads 2 unsigned bytes and orders by most significant byte first
  #
  #== Parameters:
  #offset:: how many bytes to offset for the read
  #filestream:: an opened File
  #
  #== Returns:
  #an int ordered by most significant byte first
  def readUnsignedByte_2(offset, filestream)
    #// most significant byte first
    #// |byte0|byte1| <= |unsigned int|
    byteArray = ""
    filestream.seek(offset, IO::SEEK_SET)
    byteArray = filestream.read(2, byteArray)
    return (byteArray.getbyte(0) << 8) | byteArray.getbyte(1)
  end

  #reads 4 unsigned bytes and orders by most significant byte first
  #
  #== Parameters:
  #offset:: how many bytes to offset for the read
  #filestream:: an opened File
  #
  #== Returns:
  #an int ordered by most significant byte first
  def readUnsignedByte_4(offset, filestream)
    byteArray = ""
    filestream.seek(offset, IO::SEEK_SET)
    byteArray = filestream.read(4, byteArray)
    #// most significant byte first
    #// |byte0|byte1|byte2|byte3| <= |unsigned int|
    return (byteArray.getbyte(0)<<24) | (byteArray.getbyte(1)<<16) | (byteArray.getbyte(2)<<8) | byteArray.getbyte(3)
  end

  #reads the data from the directory
  #
  #== Parameters:
  #dataOffset:: how many bytes to offset
  #numElements:: number of elements in the file computed by gatherInformation
  #filestream:: an opened File
  #
  #== Returns:
  #an array of arrays, each with information from the directory
  #[name, tag number, element type, element size, number of elements, data size, data offset]
  def readDirectoryEntry(filestream, dataOffset, numElements)
    filestream.seek(dataOffset, IO::SEEK_SET)
    byteArray = ""
    filestream.read(28*numElements, byteArray)
    directory = []
    pos = -1


    #directory structure
    #[name, tag number, element type, element size, number of elements, data size, data offset]
    (0..(numElements-1)).each do |i|
      directory[i] = []
      #// name
      name = ""
      name << byteArray.getbyte(pos+=1).chr
      name << byteArray.getbyte(pos+=1).chr
      name << byteArray.getbyte(pos+=1).chr
      name << byteArray.getbyte(pos+=1).chr
      directory[i] << name
      #// tag number
      tag_number = byteArray.getbyte(pos+=1)<<24 | byteArray.getbyte(pos+=1)<<16 | byteArray.getbyte(pos+=1)<<8 | byteArray.getbyte(pos+=1)
      directory[i] << tag_number
      #// element type
      element_type = byteArray.getbyte(pos+=1)<<8 | byteArray.getbyte(pos+=1)
      directory[i] << element_type
      #// element size
      element_size = byteArray.getbyte(pos+=1)<<8 | byteArray.getbyte(pos+=1)
      directory[i] << element_size
      #// number of elements
      number_of_elements = byteArray.getbyte(pos+=1)<<24 | byteArray.getbyte(pos+=1)<<16 | byteArray.getbyte(pos+=1)<<8 | byteArray.getbyte(pos+=1)
      directory[i] << number_of_elements
      #// data size
      data_size = byteArray.getbyte(pos+=1)<<24 | byteArray.getbyte(pos+=1)<<16 | byteArray.getbyte(pos+=1)<<8 | byteArray.getbyte(pos+=1)
      directory[i] << data_size
      #// data offset
      data_offset = byteArray.getbyte(pos+=1)<<24 | byteArray.getbyte(pos+=1)<<16 | byteArray.getbyte(pos+=1)<<8 | byteArray.getbyte(pos+=1)
      directory[i] << data_offset
      #// we do not save the dataHandle field
      pos += 4;
    end
    return directory
  end


  #directory structure
  #[name, tag number, element type, element size, number of elements, data size, data offset]
  #this is for easier index into the each directory array
  #
  #== Parameters:
  #array:: an array with information from the directory
  #element:: a string with type of information from the directory to retrieve: [name, tag_number, element_type, element_size, number_of_elements, data_size, data_offset
  #
  #== Returns:
  #the element from the array
  def get(array, element)
    if element == "name"
      return array[0]
    elsif element == "tag_number"
      return array[1]
    elsif element == "element_type"
      return array[2]
    elsif element == "element_size"
      return array[3]
    elsif element == "number_of_elements"
      return array[4]
    elsif element == "data_size"
      return array[5]
    elsif element == "data_offset"
      return array[6]
    else
      return array[0]
    end
  end

  #counts the number of samples and number of bases contained in this ABIF file
  #
  #== Parameters:
  #directory:: an array of array generated from readDirectoryEntry
  #numElements:: an int indicating the number of elements in this ABIF file
  #
  #== Returns:
  #number of samples and number of bases contained in this ABIF file
  def gatherInformation(directory, numElements)
    numSamples = 0
    numBases = 0

    (0..(numElements-1)).each do |i|
      if (get(directory[i],"name") == "DATA") && (get(directory[i], "tag_number") == 9)
        numSamples = get(directory[i], "number_of_elements") #number of elements
      else
        if (get(directory[i], "name") == "PBAS") && (get(directory[i], "tag_number") == 2)
          numBases = get(directory[i], "number_of_elements") #number of elements
        end
      end
    end

    return numSamples, numBases
  end

  #extracts the trace information for the bases
  #
  #== Parameters:
  #filestream:: an open File
  #directory:: an array of array generated by readDirectoryEntry
  #numElements:: an int indicating the number of elements in this ABIF file
  #numSamples:: an int calculated by gatherInformation
  #
  #== Returns:
  #four arrays with trace data in the order ACGT
  def getSamples(filestream, directory, numElements, numSamples)
    samples_a = []
    samples_c = []
    samples_g = []
    samples_t = []

    #// we guess the order being GATC, as Ferreira and Staden does
    (0..numElements-1).each do |i|
      tag_number = get(directory[i], "tag_number")
      if (get(directory[i],"name") == "DATA") && ([9,10,11,12].include? tag_number)
        byteArray_samples = ""
        filestream.seek(get(directory[i],"data_offset"), IO::SEEK_SET)
        filestream.read(get(directory[i], "number_of_elements")*2, byteArray_samples)
        pos = -1
        if tag_number == 9 #G
          (0..numSamples-1).each do |j|
            value = byteArray_samples.getbyte(pos+=1) << 8 | byteArray_samples.getbyte(pos+=1)
            samples_g[j] = value
          end
        elsif tag_number == 10 #A
          (0..numSamples-1).each do |j|
            value = byteArray_samples.getbyte(pos+=1) << 8 | byteArray_samples.getbyte(pos+=1)
            samples_a[j] = value
          end
        elsif tag_number == 11 #T
          (0..numSamples-1).each do |j|
            value = byteArray_samples.getbyte(pos+=1) << 8 | byteArray_samples.getbyte(pos+=1)
            samples_t[j] = value
          end
        else #C
          (0..numSamples-1).each do |j|
            value = byteArray_samples.getbyte(pos+=1) << 8 | byteArray_samples.getbyte(pos+=1)
            samples_c[j] = value
          end
        end
      end
    end
    return samples_a, samples_c, samples_g, samples_t
  end

  #extracts the called sequence information
  #
  #== Parameters:
  #filestream:: an open File
  #directory:: an array of array generated by readDirectoryEntry
  #numElements:: an int indicating the number of elements in this ABIF file
  #numBases:: an int calculated by gatherInformation
  #
  #== Returns:
  #an array with the called sequence
  def getCalledSequence(filestream, directory, numElements, numBases)
    calledSequence = []
    (0..numElements-1).each do |i|
      if (get(directory[i], "name") == "PBAS") && (get(directory[i], "tag_number") == 2)
        byteArray_seq = ""
        filestream.seek(get(directory[i], "data_offset"))
        filestream.read(numBases,byteArray_seq)
        (0..numBases-1).each do |j|
          calledSequence[j] = byteArray_seq.getbyte(j).chr
        end
      end
    end
    return calledSequence
  end

  #extracts the quality score associated with the called sequence
  #
  #== Parameters:
  #filestream:: an open File
  #directory:: an array of array generated by readDirectoryEntry
  #numElements:: an int indicating the number of elements in this ABIF file
  #numBases:: an int calculated by gatherInformation
  #
  #== Returns:
  #an array with the quality scores
  def getQualityScores(filestream, directory, numElements, numBases)
    qualityScore = []
    (0..numElements-1).each do |i|
      if (get(directory[i], "name") == "PCON") && (get(directory[i], "tag_number") == 2)
        byteArray_seq = ""
        filestream.seek(get(directory[i], "data_offset"))
        filestream.read(numBases,byteArray_seq)
        (0..numBases-1).each do |j|
          qualityScore[j] = byteArray_seq.getbyte(j)
        end
      end
    end
    return qualityScore
  end

  #extracts the trace information for the bases
  #
  #== Parameters:
  #filestream:: an open File
  #directory:: an array of array generated by readDirectoryEntry
  #numElements:: an int indicating the number of elements in this ABIF file
  #numBases:: an int calculated by gatherInformation
  #
  #== Returns:
  #an array with the indexes of the peaks
  def getPeakIndexes(filestream, directory, numElements, numBases)
    peakIndexes = []
    (0..numElements-1).each do |i|
      if (get(directory[i], "name") == "PLOC") && (get(directory[i], "tag_number") == 2)
        byteArray_peak = ""
        filestream.seek(get(directory[i], "data_offset"), IO::SEEK_SET)
        filestream.read(get(directory[i], "number_of_elements")*4, byteArray_peak)
        pos = -1
        (0..numBases-1).each do |j|
          peakIndex = byteArray_peak.getbyte(pos+=1) << 8 | byteArray_peak.getbyte(pos+=1)
          peakIndexes[j] = peakIndex
        end
      end
    end
    return peakIndexes
  end

end
