#ifndef __MESH3D__LIST_H__
#define __MESH3D__LIST_H__

template<class T>
struct Node {
	T data;
	Node *next;
};

template<class T>
Node<T> *list_allocate(int size) {
	return new Node<T>[size];
}

template<class T>
void list_deallocate(Node<T> *p) {
	delete[] p;
}

template<class T>
Node<T> *addNode(const T &data, Node<T> **mem = 0) {
	if (mem == 0) {
		Node<T> *p = new Node<T>();
		mem = &p;
	}
	(*mem)->data = data;
	(*mem)++;
	return ((*mem)-1);
}

template<class T>
void addHead(const T &data, Node<T> **head, Node<T> **mem = 0) {
	Node<T> *q = addNode(data, mem);
	q->next = *head;
	*head = q;
}

#endif
