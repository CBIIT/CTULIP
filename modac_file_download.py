import utils
import os

model_url = "https://modac.cancer.gov/api/v2/dataObject/NCI_DOE_Archive/NCI-DOE_RT/Canine_tumor_classifier/Tumor_type_classifier_models"

cnn_17 = "17CT.model.h5"
cnn_18 = "18CT.model.h5"

models = [cnn_17, cnn_18]

for model in models:
    data_loc = utils.fetch_file(model_url + model, unpack=False, md5_hash=None, subdir="models")
    print('Data downloaded and stored at: ' + data_loc)
    data_path = os.path.dirname(data_loc)
    print(data_path)
