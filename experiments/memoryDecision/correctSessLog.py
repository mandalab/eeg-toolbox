# Correct memoryDecision session logs for incorrect timestamps for STUDY_FEEDBACK_START events

import os

sub_sess = {'NIH050': [0,1,2,3,4]}

for subj, sessions in sub_sess.items():
    for i in sessions:
    
        print i    
        # first rename the original session log
        vanilla = "/Users/zaworaca/Documents/patient_extract/{}/behavioral/memoryDecision/session_{}/session.log".format(subj,sessions[i])
        orig = "/Users/zaworaca/Documents/patient_extract/{}/behavioral/memoryDecision/session_{}/session_orig.log".format(subj,sessions[i])

        print(vanilla)
        print(orig)

        os.rename(vanilla,orig)

        fout = open(vanilla,'w')
    
        with open(orig, "r") as f:
            #content = f.readlines()
            for line in f:
                values = line.split("\t")
                #print len(values)
                if len(values)>5:
                    if values[4]=='STUDY_FEEDBACK_START':
                        currentTS = values[0]
                        #print len(currentTS)
                        # Check if this is an incorrect timestamp
                        if len(currentTS) < 13:
                            lastVal = lastLine.split("\t")
                            correctedTS = str(int(lastVal[0]))
                            values[0] = correctedTS
                # Write corrected values into new file
                fout.write('\t'.join(map(str,values)))
                lastLine = line
        
        fout.close()