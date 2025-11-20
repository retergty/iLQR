#pragma once
#include <iterator>

template <typename T>
class IntrusiveList;

template <typename Derive>
class IntrusiveListNode
{
public:
  IntrusiveListNode() = default;

  // TO DO: if it is in List, remove it
  virtual ~IntrusiveListNode() = default;

private:
  friend class IntrusiveList<Derive>;

  IntrusiveListNode *prev_{nullptr};
  IntrusiveListNode *next_{nullptr};
  IntrusiveList<Derive> *owner_{nullptr};
};

template <typename T>
class IntrusiveList
{
public:
  IntrusiveList()
  {
    end_.owner_ = this;
    end_.next_ = &end_;
    end_.prev_ = &end_;
    head_ = &end_;
  }
  // TO DO: clear all node ownerness
  ~IntrusiveList() = default;

  class iterator
  {
  public:
    explicit iterator(IntrusiveListNode<T> *now = nullptr) : now_(now) {};
    operator bool() const
    {
      return now_ == nullptr;
    }
    T &operator*()
    {
      return static_cast<T &>(*now_);
    }
    T *operator->()
    {
      return static_cast<T *>(now_);
    }

    iterator &operator++()
    {
      now_ = now_->next_;
      return *this;
    }

    iterator operator++(int)
    {
      iterator tmp = *this;
      now_ = now_->next_;
      return tmp;
    }
    iterator &operator--()
    {
      now_ = now_->prev_;
      return *this;
    }
    iterator operator--(int)
    {
      iterator tmp = *this;
      tmp.now_ = now_->prev_;
      return *this;
    }
    bool operator==(const iterator &rhs)
    {
      return now_ == rhs.now_;
    }
    bool operator!=(const iterator &rhs)
    {
      return now_ != rhs.now_;
    }
    iterator operator+(int index)
    {
      iterator tmp = *this;
      for (int i = 0; i < index; ++i)
      {
        tmp = tmp->next_;
      }
      return tmp;
    }
    iterator &operator+=(int index)
    {
      for (int i = 0; i < index; ++i)
      {
        now_ = now_->next_;
      }
      return *this;
    }

    friend class IntrusiveList<T>;

  private:
    IntrusiveListNode<T> *now_;
  };

  /**
   * Remove the node from List
   * @param [in] Node: point to node
   * @return -1: error, 0 success
   */
  int remove(IntrusiveListNode<T> const *node_ptr)
  {
    IntrusiveList *owner = node_ptr->owner_;
    if (owner == nullptr)
    {
      return -1;
    }

    IntrusiveListNode<T> *last = node_ptr->prev_;
    IntrusiveListNode<T> *next = node_ptr->next_;

    last->next_ = next;
    next->prev_ = last;

    // if it is head, change it
    if (head_ == node_ptr)
    {
      head_ = next;
    }

    node_ptr->owner_ = nullptr;

    return 0;
  }

  /**
   * insert the node to this List ahead pos
   * @param [in] pos: insert position
   * @param [in] Node: node to be inserted
   * @return iterator point to inserted node
   */
  iterator insert(iterator &pos, IntrusiveListNode<T> &node)
  {
    if (!node.owner_)
    {
      // this node is in other list, remove it first
      remove(&node);
    }

    IntrusiveListNode<T> *insert_ptr = pos.now_;

    node.next_ = insert_ptr;
    node.prev_ = insert_ptr->prev_;
    node.owner_ = this;

    insert_ptr->prev_->next_ = &node;
    insert_ptr->prev_ = &node;

    if (insert_ptr == head_)
    {
      head_ = &node;
    }

    iterator ret{&node};
    return ret;
  }

  /**
   * Erases the specified elements from the container.
   * @param [in] pos: removes the element at pos
   * @return Iterator following the last removed element.If pos refers to the last element, then the end() iterator is returned.
   */
  iterator erase(iterator &pos)
  {
    IntrusiveListNode<T> *node_ptr = pos.now_;

    IntrusiveListNode<T> *last = node_ptr->prev_;
    IntrusiveListNode<T> *next = node_ptr->next_;

    last->next_ = next;
    next->prev_ = last;

    // if it is head, change it
    if (head_ == node_ptr)
    {
      head_ = next;
    }

    node_ptr->owner_ = nullptr;

    return iterator(next);
  }

  iterator begin() const
  {
    return iterator(head_);
  }

  iterator end() const
  {
    return iterator(const_cast<IntrusiveListNode<T> *>(&end_));
  }

private:
  IntrusiveListNode<T> *head_;
  IntrusiveListNode<T> end_;
};