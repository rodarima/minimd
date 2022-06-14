void
chan_send(Chan *chan, enum pb_req req)
{
    packbuf_send(&chan->pb[PB_SEND], req);
}

void
chan_recv(Chan *chan, enum pb_req req)
{
    packbuf_recv(&chan->pb[PB_RECV], req);
}

void
chan_wait1(Chan *chan, enum pb_req req)
{
    PackBuf *pb = &chan->pb[PB_RECV];
    packbuf_waitn(&pb, 1, req);
}
