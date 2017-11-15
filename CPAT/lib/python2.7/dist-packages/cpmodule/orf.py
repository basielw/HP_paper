import sys

class ORFFinder:
  """Find the longest ORF in a given sequence 
  
   "seq" is a string, if "start" is not provided any coden can be the start of 
   and ORF. If muliple ORFs have the longest length the first one encountered
   is printed 
   """
  def __init__(self, seq):
    self.seq = seq.upper()
    self.result = ("+",0,0,0,0)
    self.winner = 0
  
  def _reverse_comp(self):
    """ Returns the complementary strand, when using the negative strand
    Only used in BED files
    """
    swap = {"A":"T", "T":"A", "C":"G", "G":"C","N":"N","X":"X"}
    return "".join(swap[b] for b in self.seq)[::-1]

  """ORIGINAL CODE"""
  
  def codens(self, frame):
    """ A generator that yields DNA in one coden blocks 
    
    "frame" counts for 0. This function yields a tuple (triplet, index) with 
    index relative to the original DNA sequence 
    """
    start = frame
    while start + 3 <= len(self.seq):
      yield (self.seq[start:start+3], start)
      start += 3 

  def run_one(self, frame_number, direction,start_coden, stop_coden):
    """ Search in one reading frame """
    coden_gen = self.codens(frame_number)  
    start_codens = start_coden
    stop_codens = stop_coden   
    while True:
      try: 
        c , index = coden_gen.next()
      except StopIteration:
        break 
      # Lots of conditions here: checks if we care about looking for start 
      # coden then that coden is not a stop
      if c in start_codens or not start_codens and c not in stop_codens:
        orf_start = index  # we'll return the result as 0-indexed
        end = False
        while True:
          try: 
            c, index = coden_gen.next()
          except StopIteration:
            end = True
          if c in stop_codens:
            end = True
          if end:
            orf_end = index + 3 # because index is relative to start of coden
            L = (orf_end - orf_start)
            if L > self.winner:
              self.winner = L
              self.result = (direction, frame_number+1, orf_start, orf_end, L)
            if L == self.winner and orf_start < self.result[2]: #if ORFs have same length, return the one that if upstream
              #self.winner = L
              self.result = (direction, frame_number+1, orf_start, orf_end, L)
              
            break
    
  def longest_orf(self,direction,start_coden=['ATG'], stop_coden=['TAG','TAA','TGA']):
    if direction == "+":
      for frame in range(3):
        self.run_one(frame, direction,start_coden, stop_coden)
      return (self.result[4], self.result[1],self.seq[self.result[2]:self.result[3]]) #CDS length, coding frame, CDS sequence

    if direction == "-":
      self.seq = self._reverse_comp()
      for frame in range(3):
        self.run_one(frame, direction,start_coden, stop_coden)
      return (self.result[4], self.result[1],self.seq[self.result[2]:self.result[3]]) #CDS length, coding frame, CDS sequence

  """NEW CODE"""

  def make_codens_list(self,frame):
    tmp_seq = self.seq[frame:]
    count = 0
    codens_list = []
    while count < len(tmp_seq):
      codens_list.append(tmp_seq[count:count+3])
      count += 3
    if len(codens_list) > 0:
      if len(codens_list[-1]) < 3:
        del codens_list[-1]
    return codens_list

  def run_all(self, frame_number, start_coden, stop_coden, min_length):
    codens_list = self.make_codens_list(frame_number)
    all_orfs = []
    for a in range(0,len(codens_list)):
      if codens_list[a] in start_coden:
        for b in range(a,len(codens_list)):
          if (codens_list[b] in stop_coden) or ((b + 1) == len(codens_list)):
            orf_seq = "".join(codens_list[a:b+1])
            if len(orf_seq) >= min_length:
              start_pos = ((a * 3) + 1 + frame_number)
              all_orfs.append( (len(orf_seq), frame_number+1, orf_seq, start_pos) )
            break
    return all_orfs
  
  def all_orfs(self, direction, start_coden=['ATG'], stop_coden=['TAG','TAA','TGA'], min_length=25, kind="ALL"):
    result = []

    # if the user only wants the longest ORF, we have to set the minimum length to 0
    if kind == "LONGEST":
      min_length = 0

    if direction == '+':
      for frame in range(3):
        for i in self.run_all(frame, start_coden, stop_coden, min_length):
          result.append(i)
    
    if direction == '-':
      self.seq = self._reverse_comp()
      for frame in range(3):
        for i in self.run_all(frame, start_coden, stop_coden, min_length):
          result.append(i)

    # if the user only wants the longest ORF, we can find it using this code
    if (kind == "LONGEST") and (result != []):
      #find the longest ORFs and put them in a list max_indices that contains the indexes of the longest ORFs in result
      start_list = [item[0] for item in result]
      maximum = max(start_list)
      max_indices = [i for i, j in enumerate(start_list) if j == maximum]
      #if ORFs have the same length, return the one that is upstream
      if len(max_indices) == 1:
        result = [ result[max_indices[0]] ]
      else:
        result = [result[i] for i in max_indices]
        result.sort(key=lambda x: x[3])
        result = [result[0]]
    
    return result  #CDS length, coding frame, CDS sequence, CDS starting position

#===================
def little_test():
  seq=''
  for line in open(sys.argv[1],'r'):
    line=line.rstrip('\n\r')
    if line.startswith('>'):
  	  continue
    seq	+= line
  (l,f,s) = ORFFinder(seq).longest_orf(sys.argv[2])
  print str(l) + '\t' + str(f) + '\t' + s
  
if __name__ == "__main__":
  little_test()