docker build --platform linux/amd64 -t my_container_image .
docker run -p 80:80 --platform linux/amd64 my_container_image
