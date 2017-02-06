import tensorflow as tf
from tensorflow.examples.tutorials.mnist import input_data

mnist=input_data.read_data_sets('./mnist_data/',one_hot=True)

batch_x, batch_y=mnist.train.next_batch(100)
print(type(batch_x),type(batch_y))