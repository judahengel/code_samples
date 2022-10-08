
from typing import Pattern, Union, Tuple, List, Dict, Any
import numpy as np
import numpy.typing as npt

from typing import Pattern, Union, Tuple, List, Dict, Any
import numpy as np
import numpy.typing as npt

"""
Some type annotations
"""
from typing import Pattern, Union, Tuple, List, Dict, Any

import numpy as np
import numpy.typing as npt

"""
Some type annotations
"""
Numeric = Union[float, int, np.number, None]

"""
Global list of parts of speech
"""
POS = ['ADJ', 'ADP', 'ADV', 'AUX', 'CCONJ', 'DET', 'INTJ', 'NOUN', 'NUM',
       'PART', 'PRON', 'PROPN', 'PUNCT', 'SCONJ', 'SYM', 'VERB', 'X']

"""
Utility functions for reading files and sentences
"""
def read_sentence(f):
    sentence = []
    while True:
        line = f.readline()
        if not line or line == '\n':
            return sentence
        line = line.strip()
        word, tag = line.split("\t", 1)
        sentence.append((word, tag))

def read_corpus(file):
    f = open(file, 'r', encoding='utf-8')
    sentences = []
    while True:
        sentence = read_sentence(f)
        if sentence == []:
            return sentences
        sentences.append(sentence)


"""
3.1: Supervised learning
Param: data is a list of sentences, each of which is a list of (word, POS) tuples
Return: P(X0), 1D array; Tprob, 2D array; Oprob, dictionary {word:probabilities} 
"""
def learn_model(data:List[List[Tuple[str]]]
                ) -> Tuple[npt.NDArray, npt.NDArray, Dict[str,npt.NDArray]]:
    POS_prob = np.zeros(17)
    transition_prob = np.zeros(shape = (17, 17))
    obs_dict = {}
    obs_dict_totals = np.zeros(17)
    counter = 0
    for sentence in data:
        first_word_POS = sentence[0][1]
        idx = POS.index(first_word_POS)
        POS_prob[idx] += 1
        for ii in range(len(sentence)-1):
            if sentence[ii][0] not in obs_dict.keys():
                obs_dict[sentence[ii][0]] = np.zeros(17)
            idx = POS.index(sentence[ii][1])
            obs_dict[sentence[ii][0]][idx] += 1
            obs_dict_totals[idx] +=1
            idx_prev = POS.index(sentence[ii][1])
            idx_next = POS.index(sentence[ii+1][1])
            transition_prob[idx_next, idx_prev] += 1
            counter+=1
        counter+=1
        if sentence[-1][0] not in obs_dict.keys():
            obs_dict[sentence[-1][0]] = np.zeros(17)
        idx = POS.index(sentence[-1][1])
        obs_dict[sentence[-1][0]][idx] += 1
        obs_dict_totals[idx] += 1
    POS_prob /= len(data)
    for key in obs_dict.keys():
        obs_dict[key] = np.divide(obs_dict[key], obs_dict_totals, out=np.zeros_like(obs_dict[key]), where=obs_dict_totals!=0)
    for jj in range(17):
        transition_prob[:, jj] = np.divide(transition_prob[:, jj],  sum(transition_prob[:, jj]), out=np.zeros_like(transition_prob[:, jj]), where=transition_prob[:, jj] != 0)
    return POS_prob, transition_prob, obs_dict



"""
3.2: Viterbi forward algorithm
Param: P(X0), 1D array; Tprob, 2D array; Oprob, dictionary {word:probabilities}; obs, list of words (strings)
Return: m, 1D array; pointers, 2D array
"""
def viterbi_forward(X0:npt.NDArray,
                    Tprob:npt.NDArray,
                    Oprob:Dict[str,npt.NDArray],
                    obs:List[str]
                    ) -> Tuple[npt.NDArray, npt.NDArray]:
  m = X0.copy()
  pointers = np.zeros(shape = (len(obs),17))
  m = X0.copy()
  m2 = X0.copy()
  for counter, observation in enumerate(obs):
      for row_state, row in enumerate(Tprob):
        max_val = 0
        max_prev_state = 0
        for prev_state in range(len(row)):
          if (row[prev_state]*m[prev_state]) > max_val:
            max_val = row[prev_state]*m[prev_state]
            max_prev_state = prev_state
        pointers[counter, row_state] = max_prev_state
        m2[row_state] = max_val
      if observation in Oprob.keys():
        m2 *= Oprob[observation]
      m = m2.copy()
  return m, pointers

"""
3.2: Viterbi backward algorithm
Param: m, 1D array; pointers, 2D array
Return: List of most likely POS (strings)
"""
def viterbi_backward(m:npt.NDArray,
                     pointers:npt.NDArray
                     ) -> List[str]:
  most_likely_path = []
  idx = np.where(m == max(m))[0][0]
  for counter in reversed(range(len(pointers))):
    most_likely_path.append(idx)
    idx = (int) (pointers[counter][idx])
  most_likely_path.reverse()
  path = [POS[i] for i in most_likely_path]
  return path


"""
3.3: Evaluate Viterbi by predicting on data set and returning accuracy rate
Param: P(X0), 1D array; Tprob, 2D array; Oprob, dictionary {word:probabilities}; data, list of lists of (word,POS) pairs
Return: Prediction accuracy rate
"""
def evaluate_viterbi(X0:npt.NDArray,
                     Tprob:npt.NDArray,
                     Oprob:Dict[str,npt.NDArray],
                     data:List[List[Tuple[str]]]
                     ) -> float:
  total_words = 0;
  correct_predictions = 0;
  for sentence in data:
    m, pointers = viterbi_forward(X0, Tprob, Oprob, [word[0] for word in sentence])
    result = viterbi_backward(m, pointers)
    for counter, word in enumerate(sentence):
      if result[counter] == word[1]:
        correct_predictions +=1
      total_words+=1

  return correct_predictions/total_words


"""
3.4: Forward algorithm
Param: P(X0), 1D array; Tprob, 2D array; Oprob, dictionary {word:probabilities}; obs, list of words (strings)
Return: P(XT, e_1:T)
"""
def forward(X0:npt.NDArray,
            Tprob:npt.NDArray,
            Oprob:Dict[str,npt.NDArray],
            obs:List[str]
            ) -> npt.NDArray:
  at = X0.copy()
  for counter, observation in enumerate(obs):
    np.matmul(Tprob, at, out=at)
    if observation in Oprob.keys():
      at *= Oprob[observation]
  return at

"""
3.4: Backward algorithm
Param: Tprob, 2D array; Oprob, dictionary {word:probabilities}; obs, list of words (strings); k, timestep
Return: P(e_k+1:T | Xk)
"""
def backward(Tprob:npt.NDArray,
             Oprob:Dict[str,npt.NDArray],
             obs:List[str],
             k:int
             ) -> npt.NDArray:
  b = np.ones(len(Tprob))
  obs = obs[k+1:]
  obs.reverse()
  for observation in obs:
    if observation in Oprob.keys():
      b *= Oprob[observation]
    b = np.matmul(Tprob.transpose(), b, out = b)
  return b

"""
3.4: Forward-backward algorithm
Param: P(X0), 1D array; Tprob, 2D array; Oprob, dictionary {word:probabilities}; obs, list of words (strings); k, timestep
Return: P(Xk | e_1:T)
"""
def forward_backward(X0:npt.NDArray,
                     Tprob:npt.NDArray,
                     Oprob:Dict[str,npt.NDArray],
                     obs:List[str],
                     k:int
                     ) -> npt.NDArray:
    forward_result = forward(X0, Tprob, Oprob, obs[0:k+1])
    backward_result = backward(Tprob, Oprob, obs, k);
    forward_backward_result = forward_result * backward_result
    forward_backward_result = forward_backward_result/sum(forward_backward_result)
    return forward_backward_result


"""
3.5: Expected observation probabilities given data sequence
Param: P(X0), 1D array; Tprob, 2D array; Oprob, dictionary {word:probabilities}; data, list of lists of words
Return: New Oprob, dictionary {word:probabilities}
"""
def expected_emissions(X0:npt.NDArray,
                       Tprob:npt.NDArray,
                       Oprob:Dict[str,npt.NDArray],
                       data:List[List[str]]
                       ) -> Dict[str,npt.NDArray]:

  expected_emissions_dict = {}
  expected_emissions_totals = np.zeros(len(X0))
  for sentence in data:
    for counter, word in enumerate(sentence):
      for_back_return = forward_backward(X0, Tprob, Oprob, sentence, counter)
      if word in expected_emissions_dict.keys():
        expected_emissions_dict[word] = (expected_emissions_dict[word] + for_back_return)
      else:
        expected_emissions_dict[word] = for_back_return
      expected_emissions_totals += for_back_return
  for w in expected_emissions_dict.keys():
    np.divide(expected_emissions_dict[w],  expected_emissions_totals, out=expected_emissions_dict[w], where=expected_emissions_totals != 0)

  return expected_emissions_dict

if __name__ == "__main__":
    # Run below for 3.3
    train = read_corpus('train.upos.tsv')
    test = read_corpus('test.upos.tsv')
    X0, T, Obs_dict = learn_model(train)
    print("Train accuracy:", evaluate_viterbi(X0, T, Obs_dict, train))
    print("Test accuracy:", evaluate_viterbi(X0, T, Obs_dict, test))
    #Run below for 3.5
    obs = [[pair[0] for pair in sentence] for sentence in [test[0]]]
    Onew = expected_emissions(X0, T, Obs_dict, obs)
    print(Onew)