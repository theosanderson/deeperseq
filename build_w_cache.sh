export DOCKER_BUILDKIT=1
    
# Use cache from remote repository, tag as latest, keep cache metadata
docker build  --platform linux/amd64  -t theosanderson/deeperseq:latest \
      --cache-from theosanderson/deeperseq:latest \
      --build-arg BUILDKIT_INLINE_CACHE=1 .
