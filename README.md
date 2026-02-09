# Overview  
This is a neural network I made from scratch in **Java**.

I used the following tutorial to learn the math and algorithms that go into building a neural network 
[Neural Network Tutorial](https://blog.stackademic.com/learn-to-build-a-neural-network-from-scratch-yes-really-cac4ca457efc)

**ImagePrep.java** converts an image into an array of 4096 grayscale values that are then used to discipher whether the image is a burger or not.
**Neural.java** houses the backpropogation and feed forward algorithms, as well as all the matrix operations written from scratch needed to create the neural network.

The neural network has around a 62% success rate at determining whether an image is or isn't a burger.

## Running this program

In order to properly run this program, ensure the following:

1) All burger, nonburger, and testburger images are properly downloaded to your computer
2) All burger images used to train and test the neural network must be named with the naming convention in the images folder of this repository
3) The path to these images must be updated on lines 314, 318, and 333 to your personal path to the image



