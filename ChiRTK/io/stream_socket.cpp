#include"gnss_core_types.h"

#pragma comment(lib,"WS2_32.lib")
#pragma warning(disable:4996)


namespace
{
    /*!
    @brief 关闭套接字句柄并释放 Winsock 资源
    @param sock 待关闭的套接字句柄引用，关闭后置为 INVALID_SOCKET
    @return 无
    @note 连接失败的各分支统一走这里，避免重复写 closesocket + WSACleanup
    */
    void ReleaseSocket(SOCKET& sock)
    {
        if (sock != INVALID_SOCKET)
        {
            closesocket(sock);
            sock = INVALID_SOCKET;
        }
        WSACleanup();
    }
}


/*!
@brief 创建 TCP 客户端套接字并连接到指定 IP 地址与端口
@param sock 输出的已连接套接字句柄引用
@param IP 目标服务器 IP 地址字符串
@param Port 目标服务器端口号
@return true=连接成功 false=套接字初始化、创建或连接失败
@note 连接成功后设置 5000 ms 接收超时；失败时关闭套接字并清理 Winsock
*/
bool ConnectStreamSocket(SOCKET& sock, const char IP[], const unsigned short Port)
{
    WSADATA wsaData{};
    if (WSAStartup(MAKEWORD(2, 2), &wsaData) != 0) return false;

    sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock == INVALID_SOCKET)
    {
        WSACleanup();
        return false;
    }

    SOCKADDR_IN address{};
    address.sin_family = AF_INET;
    address.sin_addr.S_un.S_addr = inet_addr(IP);
    address.sin_port = htons(Port);

    if (connect(sock, reinterpret_cast<SOCKADDR*>(&address), sizeof(address)) == SOCKET_ERROR)
    {
        ReleaseSocket(sock);
        return false;
    }

    // 接收超时设为 5 秒，避免阻塞在 recv 上
    const int receiveTimeoutMs = 5000;

    if (setsockopt(sock, SOL_SOCKET, SO_RCVTIMEO,
        reinterpret_cast<const char*>(&receiveTimeoutMs), sizeof(receiveTimeoutMs)) == SOCKET_ERROR)
    {
        ReleaseSocket(sock);
        return false;
    }
    return true;
}

/* Close the network connection. */
/*!
@brief 关闭已建立的套接字连接并释放 Winsock 资源
@param sock 待关闭的套接字句柄引用
@return 无
@note 关闭后套接字句柄不再可用，需重新调用 ConnectStreamSocket 才能再次通信
*/
void DisconnectStreamSocket(SOCKET& sock)
{
    ReleaseSocket(sock);
}
