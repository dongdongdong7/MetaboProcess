import NeatMS as ntms
from  collections import Counter
def neatms_predict(raw_data_folder_path, feature_table_path, model_path, res_path, threshold = 0.22):
  input_data = "xcms"
  experiment = ntms.Experiment(raw_data_folder_path, feature_table_path, input_data)
  exp = experiment
  sizes = []
  print("# Feature collection:",len(exp.feature_tables[0].feature_collection_list))
  
  for consensus_feature in exp.feature_tables[0].feature_collection_list:
    sizes.append(len(consensus_feature.feature_list))

  c = Counter(sizes)
  print("Number of consensus features:")
  for size, count in c.most_common():
    print("   of size %2d : %6d" % (size, count))
  print("        total : %6d" % len(exp.feature_tables[0].feature_collection_list)) 
  
  nn_handler = ntms.NN_handler(experiment)
  nn_handler.create_model(model=model_path)
  nn_handler.get_model_summary()
  nn_handler.predict_peaks(threshold)
  
  exp = experiment
  hq_sizes = []
  lq_sizes = []
  n_sizes = []
  sizes = []
  print("# Feature collection:",len(exp.feature_tables[0].feature_collection_list))
  for consensus_feature in exp.feature_tables[0].feature_collection_list:
      hq_size = 0
      lq_size = 0
      n_size = 0
      for feature in consensus_feature.feature_list:
          for peak in feature.peak_list:
              if peak.valid:
                  if peak.prediction.label == "High_quality":
                      hq_size += 1
                  if peak.prediction.label == "Low_quality":
                      lq_size += 1
                  if peak.prediction.label == "Noise":
                      n_size += 1

      hq_sizes.append(hq_size)
      lq_sizes.append(lq_size)
      n_sizes.append(n_size)
      sizes.append(len(consensus_feature.feature_list))
      
  c = Counter(hq_sizes)
  print("\nNumber of consensus features labeled as 'High quality':")
  for size, count in c.most_common():
      print("   of size %2d : %6d" % (size, count))
  print("        total : %6d" % len(exp.feature_tables[0].feature_collection_list))

  c = Counter(lq_sizes)
  print("\nNumber of consensus features labeled as 'Low quality':")
  for size, count in c.most_common():
      print("   of size %2d : %6d" % (size, count))
  print("        total : %6d" % len(exp.feature_tables[0].feature_collection_list))

  c = Counter(n_sizes)
  print("\nNumber of consensus features labeled as 'Noise':")
  for size, count in c.most_common():
      print("   of size %2d : %6d" % (size, count))
  print("        total : %6d" % len(exp.feature_tables[0].feature_collection_list))
  
  experiment.export_csv(res_path, export_classes = ["High_quality", "Low_quality"])
