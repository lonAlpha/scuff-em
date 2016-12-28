#include "scuff-scatter-adapt.h"

#include <string.h>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <pthread.h>
#include <errno.h>

int OpenListenFd(int port)
{
  int listenfd, optval(1);
  struct sockaddr_in serveraddr;
  if ((listenfd=socket(AF_INET, SOCK_STREAM,0)) < 0) {
    fprintf(stderr, "socket error");
    fprintf(stderr, "%d: %s\n", errno, strerror(errno));
    return -1;
  }
  if (setsockopt(listenfd, SOL_SOCKET, SO_REUSEADDR, (const void*)&optval, sizeof(int))<0) {
   return -1;
  }
  memset((char *)&serveraddr, 0, sizeof(serveraddr));
  serveraddr.sin_family = AF_INET;
  serveraddr.sin_addr.s_addr=htonl(INADDR_ANY);
  serveraddr.sin_port = htons((unsigned short)port);
  if (bind(listenfd, (struct sockaddr*)&serveraddr, sizeof(serveraddr))<0) {
    fprintf(stderr, "bind error\n");
    fprintf(stderr, "%d: %s\n", errno, strerror(errno));
    return -1;
  }
  if (listen(listenfd, 128)<0) {
    fprintf(stderr, "listen error\n");
    fprintf(stderr, "%d: %s\n", errno, strerror(errno));
    return -1;
  }

  return listenfd;
}

double GetMeshSize(const double x[3], MSData *msdata)
{
  double size;
  if (msdata->firstIteration) {
    size = msdata->meshSize;
  } else {
    switch (msdata->type) {
      case RTUniform: 
        size = msdata->meshSize;
        break;
      case RTCurrent: 
        break;
      default: 
        fprintf(stderr, "Wrong RefineType\n");
        ErrExit("Wrong RefineType");
    }
  }
  return size;
}

void *MeshSizeThread(void *arg)
{
  struct sockaddr_in clientaddr;
  unsigned int clientlen;
  int listenfd, connfd;
  double x[3], size;
  MSData *msdata;
  msdata = (MSData *)arg;
  msdata->threadId = pthread_self();
  if ((listenfd = OpenListenFd(MESHPORT))<0) {
    fprintf(stderr, "Open port %d failed\n", MESHPORT);
    ErrExit("Open port failed");
  }
  while (1) {
    clientlen = sizeof(clientaddr);
    if ((connfd = accept(listenfd, (struct sockaddr*)&clientaddr, &clientlen))<0) {
      fprintf(stderr, "Accept failed\n");
      ErrExit("Accept failed");
    }
    read(connfd, x, 3*sizeof(x[0]));
    size = GetMeshSize(x, msdata);
    write(connfd, &size, sizeof(size));
    close(connfd);
  }
}

void StartBGMeshService(MSData *msdata)
{
  msdata->firstIteration = true;
  pthread_t tid;
  pthread_create(&tid, NULL, MeshSizeThread, (void*)msdata);
}

void UpdateBGMeshService(MSData *msdata)
{
  msdata->firstIteration = false;
  switch (msdata->type) {
    case RTUniform: 
      msdata->meshSize /= (1+msdata->percentage);
      break;
    default: 
      ErrExit("Error in UpdateBGMeshService");
  }
}

void EndBGMeshService(MSData *msdata)
{
  pthread_cancel(msdata->threadId);
}


