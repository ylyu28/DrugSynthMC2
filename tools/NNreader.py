import numpy as np
import tensorflow as tf
from dataclasses import dataclass
from typing import List

@dataclass
class Prediction:
    label: List[int]
    confidence: List[float]

class Model:
    def __init__(self, bundle, input_op, input_index, output_op, output_index):
        self.bundle = bundle
        self.input_op = input_op
        self.input_index = input_index
        self.output_op = output_op
        self.output_index = output_index
    
    @classmethod
    def load_model(cls, export_dir:str):

        bundle = tf.saved_moel.load(export_dir)
        signatures = bundle.signatures[tf.saved_model.DEFAULT_SERVING_SIGNATURE_DEF_KEY]

        input_info = signatures.structured_input_signature[1]['embedding_input']
        output_info = signatures.structured_outputs['time_distributed']

        input_op = input_info.name
        output_op = output_info.name

        input_index = 0
        output_index = 0

        return cls(bundle, input_op, input_index, output_op, output_index)
    
    def predict(self, SMILES:List[str]) -> Prediction:
        val = ["\n", "&", "C", "(", ")", "1", "=", "2", "O", "N", "3", "F", "[C@@H]", "#", 
               "S", "L", "[O-]", "[C@H]", "[NH+]", "[C@]", "Br", "/", "[NH3+]", "W", "4", 
               "[NH2+]", "[C@@]", "[N+]", "\\", "M", "[S@]", "5", "[N-]", "[S@@]", "[S-]", 
               "6", "7", "I", "P", "[OH+]", "[NH-]", "[P@@H]", "[P@@]", "[PH2]", "[P@]", 
               "[P+]", "[S+]", "[O+]", "[CH2-]", "[CH-]", "[SH+]", "[PH+]", "[PH]", "8", 
               "[S@@+]"]
        
        get_int = [1.0]

        for c in SMILES:
            if c in val:
                get_int.append(float(val.index(c)))
        
        while len(get_int) != 81:
            get_int.append(0.0)
        
        input_tensor = tf.convert_to_tensor(np.array(get_int).reshape(1, 81), dtype=tf.float32)

        output = self.bundle.signatures[tf.saved_model.DEFAULT_SERVING_SIGNATURE_DEF_KEY](input_tensor)

        output_values = output['time_distributed'].numpy()

        confidence = []
        label = []
        max_lab = 0
        max_conf = 0.0
        for i in range(output_values.shape[2]):
            conf_value = output_values[0, len(SMILES), i]
            confidence.append(conf_value)
            if conf_value > max_conf:
                max_conf = conf_value
                max_lab = i
            
            label.append(i)
        
        return Prediction(label, confidence)
    