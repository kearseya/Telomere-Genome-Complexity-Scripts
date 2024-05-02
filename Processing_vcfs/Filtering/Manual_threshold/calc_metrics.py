import os

### calculate the precision and recall after manually classifying variants

indir = "/path/to/project/variants/processing_vcfs/manual_filtering/hawk"

for types in os.listdir(indir):
    types_path = os.path.join(indir, types)
    for p in os.listdir(types_path):
        print("Type: ", types, "\tProb: ", p)
        tp = len(os.listdir(os.path.join(types_path, p, "TP")))
        tn = len(os.listdir(os.path.join(types_path, p, "TN")))
        fp = len(os.listdir(os.path.join(types_path, p, "FP")))
        fn = len(os.listdir(os.path.join(types_path, p, "FN")))
        if tp > 0:
            print("Precision :", tp/(tp+fp))
            print("Recall    :", tp/(tp+fn))
        else:
            print("Classify SVs, missing TPs")
