docker build --platform linux/amd64 -t my_container_image .
docker tag my_container_image theosanderson/deepseqer:latest
docker push theosanderson/deepseqer:latest