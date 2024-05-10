SPLICEBERT_PATH = "/path/to/SpliceBERT/models/model_folder"  # set the path to the folder of pre-trained SpliceBERT
import torch
from transformers import AutoTokenizer, AutoModel, AutoModelForMaskedLM, AutoModelForTokenClassification

# load tokenizer
tokenizer = AutoTokenizer.from_pretrained(SPLICEBERT_PATH)

# prepare input sequence
seq = "ACGUACGuacguaCGu"  ## WARNING: this is just a demo. SpliceBERT may not work on sequences shorter than 64nt as it was trained on sequences of 64-1024nt in length
seq = ' '.join(list(seq.upper().replace("U", "T"))) # U -> T and add whitespace
input_ids = tokenizer.encode(seq) # N -> 5, A -> 6, C -> 7, G -> 8, T(U) -> 9. NOTE: a [CLS] and a [SEP] token will be added to the start and the end of seq
input_ids = torch.as_tensor(input_ids) # convert python list to Tensor
input_ids = input_ids.unsqueeze(0) # add batch dimension, shape: (batch_size, sequence_length)

# use huggerface's official API to use SpliceBERT
# get nucleotide embeddings (hidden states)
model = AutoModel.from_pretrained(SPLICEBERT_PATH) # load model
last_hidden_state = model(input_ids).last_hidden_state # get hidden states from last layer
hiddens_states = model(input_ids, output_hidden_states=True).hidden_states # hidden states from the embedding layer (nn.Embedding) and the 6 transformer encoder layers

# get nucleotide type logits in masked language modeling
model = AutoModelForMaskedLM.from_pretrained(SPLICEBERT_PATH) # load model
logits = model(input_ids).logits # shape: (batch_size, sequence_length, vocab_size)

# finetuning SpliceBERT for token classification tasks
model = AutoModelForTokenClassification.from_pretrained(SPLICEBERT_PATH, num_labels=3) # assume the class number is 3, shape: (batch_size, sequence_length, num_labels)

# finetuning SpliceBERT for sequence classification tasks
model = AutoModelForSequenceClassification.from_pretrained(SPLICEBERT_PATH, num_labels=3) # assume the class number is 3, shape: (batch_size, sequence_length, num_labels)