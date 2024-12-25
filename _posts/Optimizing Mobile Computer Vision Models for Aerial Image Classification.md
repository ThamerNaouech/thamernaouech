---
date: '2024-12-19T11:50:54.000Z'
title: Optimizing Mobile Computer Vision Models for Aerial Image Classification
tagline: 
preview: >-
  
image: >-
  https://eod-grss-ieee.com/uploads/science/1659023090_SeasoNet.jpg
---

# Introduciton

Mobile computer vision models have revolutionized how we interact with and interpret visual data, especially on resource-constrained devices. From real-time image classification to object detection, these models are indispensable in industries like healthcare, agriculture, and environmental monitoring. 
In this blog post, I’ll walk you through my journey of creating a deep learning project that fine-tunes state-of-the-art mobile vision models for aerial image classification using the EuroSat dataset, as part of my participation in the Mobile Computer Vision Lab Course @ TUM. 

# Project Overview

The objective of this project was to classify satellite imagery from the EuroSat dataset into five categories: Highway, Industrial, Sea/Lake, Forest, and River. To achieve this, I fine-tuned multiple mobile-friendly deep learning models, including MobileNetV2, MobileNetV3, and EfficientNet-Lite. These models were chosen for their balance between computational efficiency and performance, making them ideal for deployment on mobile or edge devices.

![Transfer Learning](https://miro.medium.com/v2/resize:fit:700/0*ENPuddcWMZAOj_xU.png)

# Tools and Technologies


**Framework:** PyTorch

**Models:** MobileNetV2, MobileNetV3, EfficientNet-Lite

**Dataset:** EuroSat (RGB images)

**Techniques:** Transfer learning, data augmentation, and comparative model training




# Step-by-Step Implementation
## 1. Dataset Preparation

The EuroSat dataset contains RGB satellite imagery that is pre-labeled for land use classification. I organized the dataset into five folders corresponding to the target classes. Additionally, I excluded unnecessary files, such as .ipynb_checkpoints, to ensure clean data processing.

```jsx
      from torchvision import datasets
      train_dataset = datasets.ImageFolder(
          root="EuroSet",
          transform=transforms.Compose([
              transforms.RandomHorizontalFlip(),
              transforms.RandomRotation(10),
              transforms.ToTensor(),
              transforms.Normalize(mean=[0.4914, 0.4822, 0.4465], std=[0.2023, 0.1994, 0.2010])
          ])
      )
```

## 2. Model Selection

To balance accuracy and efficiency, I chose:
 * MobileNetV2 and MobileNetV3: Lightweight architectures optimized for mobile devices.
 * EfficientNet-Lite: A streamlined version of EfficientNet designed for resource-constrained environments.
These models were initialized with pre-trained weights from ImageNet to leverage transfer learning.
## 3. Fine-Tuning the Models

The models were modified to accommodate the five output classes of the EuroSat dataset by replacing the final fully connected layers.

## 4. Training Strategy
I implemented data augmentation techniques such as random rotations and horizontal flips to increase the diversity of the training dataset. The dataset was split into training and validation sets (80/20 split), and I used the Adam optimizer with cross-entropy loss.
```jsx
      criterion = nn.CrossEntropyLoss()
      optimizer = optim.Adam(model.parameters(), lr=0.001)
```

## 5. Comparative Analysis
Each model was trained for 10 epochs, and their validation performance was compared in terms of loss and accuracy. Here’s a summary of the process:

 * Without transformations: Baseline performance of the models.

 * With transformations: Performance improvement with augmented data.

 
## 6. Evaluation and Best Model Selection
Validation accuracy and loss were plotted to visualize the training process. The best-performing model was saved for deployment.
```jsx
      torch.save(model.state_dict(), "best_model.pth")
```
# Results and Insights

 * MobileNetV3 achieved the highest validation accuracy, showcasing its effectiveness for mobile deployments.

 * Data augmentation significantly improved model robustness, particularly for classes with fewer samples.

 * EfficientNet-Lite, while slightly heavier than MobileNet variants, offered a good balance between accuracy and computational efficiency.



# Challenges and Solutions

 * Handling Imbalanced Data: Some classes in the EuroSat dataset had fewer samples. I addressed this by using data augmentation to generate more diverse examples for these classes.

 * Model Overfitting: Dropout layers and early stopping were employed to prevent overfitting during training.

# Complete code
```jsx
    TBD
```



# References
Helber, P., Bischke, B., Dengel, A., and Borth, D., “EuroSAT: A Novel Dataset and Deep Learning Benchmark for Land Use and Land Cover Classification”, *arXiv e-prints*, Art. no. arXiv:1709.00029, 2017. doi:10.48550/arXiv.1709.00029.
