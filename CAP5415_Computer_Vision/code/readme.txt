This project performs object detection using Yolov3 on images and videos.
The contributors are: Ghayoor Shah, Mahdi Razzaghpour

First, if you want to evaluate any model on image/video, you need to run the demo_image.py/ demo_video.py file
For example, if you want to evaluate the pre-trained Yolov3 model using pre-trained Yolov3 weights,
then you need to run the following commands

python demo_image.py --weights_path requirements/yolov3.weights --image data/inssbruck.png --cfg config/yolov3_default.config
OR
python demo_video.py --weights_path requirements/yolov3.weights --image data/video1.avi --cfg config/yolov3_default.config

On the other hand, if you want to evaluate the re-trained model, you can run the following commands:

python demo_image.py --ckpt checkpoints/snapshot6000.ckpt --image data/inssbruck.png --cfg config/yolov3_default_test1.config
OR
python demo_video.py --ckpt checkpoints/snapshot6000.ckpt --image data/video1 --cfg config/yolov3_default_test1.config


If you wish to re-train the model using yolov3 weights and pre-trained yolov3 configurations, run the following command:

python train.py --weights_path requirements/yolov3.weights --cfg config/yolov3_default_test1.config

On the other hand, if you wish to re-train the model using Darknet backbone as initial weights and yolov3 configurations, run the following commands:

python train.py --weights_path requirements/darknet53.conv.74 --cfg config/yolov3_default_test1.config

Note: Due to the size constraints when uploading the file on webcourses, we have not included the weights and checkpoints in this folder.
      However, if required, we can send as onedrive link.
